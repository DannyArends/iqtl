/**
 * \file regressionsupport.h
 * \brief Header file, for regressionsupport.cpp
 *
 * last modified Sep, 2010
 * first written Apr, 2010 
 * Copyright (c) 2010 Danny Arends
 * 
 **/

 #ifndef REGRESSIONSUPPORT_H_
  #define REGRESSIONSUPPORT_H_
extern "C"
{  
  #include "datastructures.h"
  
  double Lnormal(double residual, double variance);
  dvector calculateparameters(uint nvariables, uint nsamples, dmatrix xt, dvector w, dvector y,int verbose);
  dmatrix translatematrix(uint nvariables, uint nsamples, dmatrix x,int verbose);
  double calculatestatistics(uint nvariables, uint nsamples, dmatrix xt, dvector xtwy, dvector y,dvector w, dvector* fit, dvector* residual,int verbose);
  double calculateloglikelihood(uint nsamples, dvector residuals,dvector w, double variance, dvector* Fy,int verbose);

}
#endif
