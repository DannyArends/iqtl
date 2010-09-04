/**
 * \file typedefs.h
 *
 * Copyright (c) 2010 Danny Arends
 * 
 **/
 

#ifndef TYPEDEFS_H_
  #define TYPEDEFS_H_
  
extern "C"
{  
  #include <ctype.h>
  #include <cstdio>
  #include <cstdlib>
  #include <R.h>
  #include <getopt.h>
  
  using namespace std;
  
  typedef unsigned int uint;
  
  typedef double**  dmatrix;
  typedef double*   dvector;
  typedef bool*     bvector;
  typedef char**    cmatrix;
  typedef char*     cvector;
  typedef int**     imatrix;
  typedef int*      ivector;
}
#endif
