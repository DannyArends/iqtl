import std.stdio;
import std.math; 
import std.array;
import std.string;
import std.conv;

string[] getStringArrayFromR(int* r_length, int* r_sizes, char** r_values){
  int[] sizes = r_sizes[0..(*r_length)];
  string values[];
  for(int s=0;s<(*r_length);s++){
    writefln("s: %s",s);
    string tmp;
    for(int c=0;c<sizes[s];c++){
	  tmp ~= r_values[s][c];
    }
    values ~= tmp;
  }
  return values;
}