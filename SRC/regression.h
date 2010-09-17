/**
 * \file regression.h
 * \brief Header file, for regression.cpp
 *
 * last modified Sep, 2010
 * first written Apr, 2010 
 * Copyright (c) 2010 Danny Arends
 * 
 **/

 #ifndef REGRESSION_H_
  #define REGRESSION_H_
extern "C"
{  
  #include "datastructures.h"

  double Lnormal(double residual, double variance);
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
  void nulllikelihoodbyem_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* verbose,double* out);
  void modellikelihoodbyem_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* verbose,double* out);
  void likelihoodbyem_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* verbose,double* out);
  double likelihoodbyem(uint nvariables,uint nsamples, dmatrix x, dvector w, dvector y,int verbose);
  dvector calculateparameters(uint nvariables, uint nsamples, dmatrix xt, dvector w, dvector y,int verbose);
  dmatrix translatematrix(uint nvariables, uint nsamples, dmatrix x,int verbose);
  double multivariateregression(uint nvariables, uint nsamples, dmatrix x, dvector w, dvector y, dvector Fy,int verbose);
  double nullmodel(uint nvariables, uint nsamples, dmatrix x, dvector w, dvector y,int verbose);
  bool LUdecomposition(dmatrix m, int dim, ivector ndx, int *d);
  void LUsolve(dmatrix lu, int dim, ivector ndx, dvector b);
  void LUinvert(dmatrix lu, dmatrix inv, int dim, int *ndx);
}
#endif
