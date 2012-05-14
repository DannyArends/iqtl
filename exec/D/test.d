/**
 * Functions exported by Rlib
 */

module qtl.plugins.renv.test;

import std.c.stdio;
import std.conv;

void main(string[] args){ }

extern (C){
  export void connected(int* x){
    if((*x)==0){
      (*x) = 1;
    }
  }
}
