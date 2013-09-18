/**
 * \file datastructures.h
 * \brief Header file, for datastructures.cpp
 *
 * last modified Sep, 2010
 * first written Apr, 2010
 * copyright (c) 2010 Danny Arends
 *
 * C functions for the iqtl package
 * Contains:  newbvector,newdvector,newcvector,newivector
 *            printdvector,printcvector,printivector
 *            newdmatrix,newcmatrix,newimatrix
 *            printdmatrix,printcmatrix,printimatrix
 *            freematrix,freevector,mycalloc
 *
 **/
 
#ifndef _DATASTRUCTURES_H
  #define _DATASTRUCTURES_H
  #include "typedefs.h"
  extern "C"
  {
  bvector newbvector(uint dim);
  dvector newdvector(uint dim);
  ldvector newldvector(uint dim);
  cvector newcvector(uint dim);
  ivector newivector(uint dim);

  void printldvector(dvector v, uint dim);
  void printdvector(dvector v, uint dim);
  void printcvector(cvector v, uint dim);
  void printivector(ivector v, uint dim);

  dmatrix newdmatrix(uint rows, uint cols);
  cmatrix newcmatrix(uint rows, uint cols);
  imatrix newimatrix(uint rows, uint cols);

  void printdmatrix(dmatrix m, uint rows, uint cols);
  void printcmatrix(cmatrix m, uint rows, uint cols);
  void printimatrix(imatrix m, uint rows, uint cols);

  void dvectortodmatrix(int nrow, int ncol, dvector in, dmatrix* out);
  void ivectortoimatrix(int nrow, int ncol, ivector in, imatrix* out);
  
  void freematrix(void **m, uint rows);
  void freevector(void *v);
  void *mycalloc(uint num, uint size);
  }
#endif
