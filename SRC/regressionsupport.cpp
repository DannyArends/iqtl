/**
 * \file regressionsupport.cpp
 * \brief Code file containing functions
 *
 * last modified Sep, 2010
 * first written Apr, 2010
 * Copyright (c) 2010 Danny Arends
 * 
 **/
 
#include "regressionsupport.h"
#include "LUdecomposition.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <limits> 
#include <cfloat>
#include <R.h>
#include <Rmath.h>

double Lnormal(double residual, double variance){
  return dnorm(residual,0,sqrt(variance),0);
}

dvector calculateparameters(uint nvariables, uint nsamples, dmatrix xt, dvector w, dvector y, int verbose){
  int d=0;
  double xtwj;
  dmatrix XtWX = newdmatrix(nvariables, nvariables);
  dvector XtWY = newdvector(nvariables);
  ivector indx = newivector(nvariables);

  if(verbose) Rprintf("calculating XtWX and XtWY\n");
  for(uint i=0; i<nsamples; i++){
    for(uint j=0; j<nvariables; j++){
      xtwj     = xt[j][i] * w[i];
      XtWY[j] += xtwj    * y[i];
      for(uint jj=0; jj<=j; jj++){
        XtWX[j][jj] += xtwj * xt[jj][i];
      }
    }
  }
 
  LUdecomposition(XtWX, nvariables, indx, &d);
  LUsolve(XtWX, nvariables, indx, XtWY);
  
  freematrix((void**)XtWX, nvariables);
  
  return XtWY;
}

dmatrix translatematrix(uint nvariables, uint nsamples, dmatrix x, int verbose){
  dmatrix Xt   = newdmatrix(nvariables,nsamples);
  if(verbose) Rprintf("calculating Xt\n");
  for(uint i=0; i<nsamples; i++){
    for(uint j=0; j<nvariables; j++){
      Xt[j][i] = x[i][j];
    }
  }
  return Xt;
}

double calculatestatistics(uint nvariables, uint nsamples, dmatrix xt, dvector xtwy, dvector y, dvector w, dvector* fit, dvector* residual,int verbose){
  double variance= 0.0;
  for (uint i=0; i<nsamples; i++){
    (*fit)[i]      = 0.0;
    (*residual)[i] = 0.0;
    for (uint j=0; j<nvariables; j++){
      (*fit)[i]     += xt[j][i] * xtwy[j];
    }
    (*residual)[i]   = y[i]-(*fit)[i];
    variance        += w[i]*pow((*residual)[i],2.0);
  }
  variance /= nsamples;
  return variance;
}

double calculateloglikelihood(uint nsamples, dvector residual,dvector w, double variance, dvector* Fy,int verbose){
  double logL  = 0.0;
  ldvector indL    = newldvector(nsamples);

  for (uint i=0; i<nsamples; i++){
    (*Fy)[i]  = Lnormal(residual[i],variance);
    indL[i]  += w[i] * (*Fy)[i];
    logL     += log(indL[i]);
  }
  freevector((void*)indL);
  
  return logL;
}
