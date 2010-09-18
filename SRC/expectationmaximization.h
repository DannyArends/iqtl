/**
 * \file expectationmaximization.h
 * \brief Header file, for expectationmaximization.cpp
 *
 * last modified Sep, 2010
 * first written Apr, 2010 
 * Copyright (c) 2010 Danny Arends
 * 
 **/

 #ifndef EMAXIMIZATION_H_
  #define EMAXIMIZATION_H_
extern "C"
{  
  #include "datastructures.h"
  #include "regression.h"

  void nullmodellikelihood_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* verbose,double* out);
  void modellikelihoodbyem_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* verbose,double* out);
  void lodscorebyem_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* verbose,double* out);
  double likelihoodbyem(uint nvariables,uint nsamples, dmatrix x, dvector w, dvector y,int verbose);
}
#endif
