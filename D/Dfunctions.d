import std.stdio;
import std.math; 
import std.array;
import std.string;
import std.conv;
import std.file;
import std.regex;
import std.datetime;
import core.time;
import core.vararg;
import std.c.stdlib;

import std.c.windows.windows;
import core.sys.windows.dll;

import rlib;
import raux;

enum BUFFERSIZE {BUFFER_16KB = 16_384,BUFFER_2MB = 2_097_152, BUFFER_4MB = 4_194_304,BUFFER_16MB = 16_777_216, BUFFER_64MB = 67_108_864, BUFFER_256MB = 268_435_456}

extern(C){
  export void Test(int* length,int* values){
    int Dvalues[];
    for(int x=0;x<(*length);x++){
      Dvalues ~= values[x];
    }
    writeln("Test: ",Dvalues);
  }
  
  export void CharTest(int* r_length, int* r_sizes, char** r_values){
	  string values[] = getStringArrayFromR(r_length, r_sizes, r_values);
    writeln("Test: ",values);
  }
    
 /*
  * Prints out a description of a tab separated file (any size).
  *
  * @param filename to load
  * @return Number of buffers needed to read in the entire file
  */
  export void loadLineFromFile(char** r_filenames, int* r_lines, int* r_header, char** r_sep, int* r_unquote, int* r_l_filename, int* r_s_filename, int* r_l_lines, int* output){
    string[] filenames = getStringArrayFromR(r_l_filename, r_s_filename, r_filenames);
    string filename = filenames[0];
    writefln("Filename %s", filename);
    int[] lines = r_lines[0..(*r_l_lines)];
    writefln("Lines %d", lines);
    bool header = (*r_header)==1;
    writefln("Header %s", header);
    char sep = r_sep[0][0];
    writefln("Sep '%s'", sep);
    bool unquote = (*r_unquote)==1;
    writefln("unquote %s", unquote);
    uint buffersize = BUFFERSIZE.BUFFER_16KB;
    
    if(!exists(filename) || !isfile(filename)){
      writefln("File '%s' does not exist", filename);
      return;
    }
    uint filesize = cast(uint)getSize(filename);
    long linecount=0;
    ubyte[] inputbuffer = new ubyte[buffersize];
    int[] headerline;
    int headersize =0;
    auto f = new File(filename,"rb");
    auto t0 = Clock.currTime();
    long buffercount = 0;
    long tabscount = 0;
    while(f.rawRead(inputbuffer)){
      foreach(int i,byte b ; inputbuffer){
        if(header && linecount ==0){
          headerline ~= b;
          headersize++;
        }
        if(linecount >3){
          writefln("Header: %s",intToString(headerline));
           //
           output = cast(int*)R_alloc(4,int.sizeof);
           output[0] = 100;
           output[1] = 200;
           output[2] = 300;
          return;
        }
        switch(cast(char)b){
          case '\n':
            linecount++;
            break;
          case '\t':
            tabscount++;
            break;
          default: break;
        }
      }
      delete inputbuffer; inputbuffer = new ubyte[buffersize];
      buffercount++;
    }
    auto t1 = Clock.currTime();
    f.close();
    writefln("Filesize %d: %d buffers %d lines, %d tabs in %d", filesize, buffercount, linecount, tabscount, (t1-t0));
  }
}
