/**********************************************************************
 * src/D/libs/r.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module libs.r;

import std.stdio;
import std.conv;

version (Windows) {
  import libs.libload;
  import std.loader;
    
  extern(C){
    double function(double, double, double, int) dnorm;
    double function(double, double, double, int, int) qf;
    void   function(char *, ...) Rprintf;

  void LoadR(){
    HXModule lib = load_library("R");
    
    load_function(dnorm)(lib,"Rf_dnorm4");
    load_function(qf)(lib,"Rf_qf");
    load_function(Rprintf)(lib,"Rprintf");
    
    writeln("Loaded R functionality");
  }
  
  }
  
}else{
  pragma(lib, "R");
  
  extern(C){
    double dnorm(double, double, double, int);
    double qf(double, double, double, int, int);
    void   Rprintf(char *, ...);
  }
}

unittest{
  writeln("Unit test " ~ __FILE__);
}
