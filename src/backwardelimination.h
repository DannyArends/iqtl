/**
 * \file backwardelimination.h
 * \brief Header file, for backwardelimination.cpp
 *
 * last modified Sep, 2010
 * first written Apr, 2010 
 * Copyright (c) 2010 Danny Arends
 * 
 **/

 #ifndef BACKWARDELIMINATION_H_
  #define BACKWARDELIMINATION_H_
extern "C"
{  
  #include "datastructures.h"

  void testregression();
  uint lowestindex(uint dim, dvector values);
  uint highestindex(uint dim, dvector values);
  void inverseF_R(int* df1,int* df2, double* alfa, double* out);
  double inverseF(int df1, int df2, double alfa);
  uint modelsize(uint nvariables, bvector model);
  void dropterm(uint nvariables, bvector model,uint which);
  dmatrix createdesignmatrix(uint nvariables,uint nsamples, dmatrix x, bvector model);
  void copybvector(uint dim,bvector origin,bvector target);
  void backwardelimination_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* verbose,double* out);
  void backwardelimination(uint nvariables,uint nsamples, dmatrix x, dvector w, dvector y,int verbose);
}
#endif
