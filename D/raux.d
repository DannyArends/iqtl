import std.stdio;
import std.math; 
import std.array;
import std.string;
import std.conv;

pure bool searchArray(T)(T[] haystack, T needle){
  foreach(T s; haystack){
    if(s==needle){
      return true;
    }
  }
  return false;
}

string[] getStringArrayFromR(int* r_length, int* r_sizes, char** r_values){
  int[] sizes = r_sizes[0..(*r_length)];
  string values[];
  for(int s=0;s<(*r_length);s++){
    string tmp;
    for(int c=0;c<sizes[s];c++){
	  tmp ~= r_values[s][c];
    }
    values ~= tmp;
  }
  return values;
}

string intToString(int[] bits){
  char[] r;
  foreach(int b;bits){
    r ~= cast(char)(b);
  }
  //r ~= '\0';
  return to!string(r);
}
