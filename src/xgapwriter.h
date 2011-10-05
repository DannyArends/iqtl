/**
 * \file xgapwriter.h
 * \brief Header file, for xgapwriter.cpp
 *
 * last modified Sep, 2010
 * first written Apr, 2009
 * copyright (c) 2009-2010 Danny Arends
 *
 * C functions for the iqtl package
 * Contains:  
 *
 **/

#ifndef _XGAPWRITER_H
  #define _XGAPWRITER_H
  extern "C"
  {
  #include <R.h>
  #include "datastructures.h"
  
  void r_writebinaryfile(char** filename, int* namelength, int* nrows, int* ncols, int* type,int** lengths, void** data);
  
  }
#endif
