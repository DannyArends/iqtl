/**
 * \file typedefs.h
 * \brief Header file, Some types used in the iqtl package
 *
 * last modified Sep, 2010
 * first written Apr, 2009
 * copyright (c) 2009-2010 Danny Arends
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

  typedef double**      dmatrix;
  typedef double*       dvector;
  typedef long double*  ldvector;
  typedef bool*         bvector;
  typedef char**        cmatrix;
  typedef char*         cvector;
  typedef int**         imatrix;
  typedef int*          ivector;
  }
#endif
