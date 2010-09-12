/**
 * \file regression.h
 * \brief Header file, for regression.cpp
 *
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
  double betacf(double a, double b, double x);
  double gammln(double xx);
  double betai(double a, double b, double x);
  double inverseF(int df1, int df2, double alfa);
  uint modelsize(uint nvariables, bvector model);
  void dropterm(uint nvariables, bvector model,uint which);
  dmatrix createdesignmatrix(uint nvariables,uint nsamples, dmatrix x, bvector model);
  void copybvector(uint dim,bvector origin,bvector target);
  void backwardelimination(uint nvariables,uint nsamples, dmatrix x, dvector w, dvector y);
  double likelihoodbyem(uint nvariables,uint nsamples, dmatrix x, dvector w, dvector y);
  double multivariateregression(uint nvariables, uint nsamples, dmatrix x, dvector w, dvector y, dvector Fy);
  bool LUdecomposition(dmatrix m, int dim, ivector ndx, int *d);
  void LUsolve(dmatrix lu, int dim, ivector ndx, dvector b);
  void LUinvert(dmatrix lu, dmatrix inv, int dim, int *ndx);
}
#endif
