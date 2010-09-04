/**
 * \file regression.h
 * \brief Header file, for regression.cpp
 *
 * Copyright (c) 2010 Danny Arends
 * 
 **/

 #ifndef REGRESSION_H_
  #define REGRESSION_H_
  
  #include "datatypes.h"

  double Lnormal(double residual, double variance);
  void testregression();
  void backwardelimination(uint nvariables,uint nsamples, dmatrix x, dvector w, dvector y);
  double likelihoodbyem(uint nvariables,uint nsamples, dmatrix x, dvector w, dvector y);
  double multivariateregression(uint nvariables, uint nsamples, dmatrix x, dvector w, dvector y, dvector Fy);
  bool LUdecomposition(dmatrix m, int dim, ivector ndx, int *d);
  void LUsolve(dmatrix lu, int dim, ivector ndx, dvector b);
  void LUinvert(dmatrix lu, dmatrix inv, int dim, int *ndx);

#endif
