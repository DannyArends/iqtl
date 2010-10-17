/**
 * \file regression.cpp
 * \brief Code file containing functions
 *
 * last modified Sep, 2010
 * first written Apr, 2010
 * Copyright (c) 2010 Danny Arends
 * 
 **/
 
#include "regression.h"
#include "regressionsupport.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <limits> 
#include <cfloat>
#include <R.h>
#include <Rmath.h>

double multivariateregression(uint nvariables, uint nsamples, dmatrix x, dvector w, dvector y, dvector Fy, bool nullmodel, ivector nullmodellayout,int verbose){
  dmatrix Xt   = translatematrix(nvariables,nsamples,x,verbose);
  dvector XtWY = calculateparameters(nvariables,nsamples,Xt,w,y,verbose);

  if(nullmodel){
    for (uint i=1; i < nvariables; i++){
      if(nullmodellayout[(i-1)] == 1){ //SHIFTED Because the nullmodel has always 1 parameter less (The first parameter estimated mean)
        XtWY[i] = 0.0;
      }
    }
  }
  
  if(verbose){
    Rprintf("Estimated parameters:\n");
    printdvector(XtWY,nvariables);
  }
  
  dvector fit       = newdvector(nsamples);
  dvector residual  = newdvector(nsamples);
  double  variance  = calculatestatistics(nvariables, nsamples, Xt, XtWY, y, w, &fit, &residual,verbose);
  double  logLQTL   = calculateloglikelihood(nsamples, residual, w, variance, &Fy, verbose);
  
  if(verbose){
    Rprintf("Estimated response:\n");
    printdvector(fit,nsamples);

    Rprintf("Residuals:\n");
    printdvector(residual,nsamples);

    Rprintf("Estimated Fy:\n");
    printdvector(Fy,nsamples);

    Rprintf("Variance: %f\n",variance);
    Rprintf("Loglikelihood QTL: %f\n",logLQTL);
  }
  
  freematrix((void**)Xt,nvariables);
  freevector((void*)XtWY);
  freevector((void*)fit);
  freevector((void*)residual);

  return logLQTL;
}

void inverseF_R(int* df1,int* df2, double* alfa, double* out){
  (*out) = inverseF((*df1), (*df2), (*alfa));
}

double inverseF(int df1, int df2, double alfa){
  double prob=0.0, minF=0.0, maxF=100.0, halfway=50.0, absdiff=1.0;
  int count=0;
  while ((absdiff>0.001)&&(count<100)){
    count++;
    halfway= (maxF+minF)/2.0;
    prob = pbeta(df2/(df2+df1*halfway), df2/2.0, df1/2.0, 1, 0);
    if (prob<alfa) maxF= halfway;
    else minF= halfway;
    absdiff= fabs(prob-alfa);
  }
  return halfway;
}

double nullmodel(uint nvariables, uint nsamples, dmatrix x, dvector w, dvector y,ivector nullmodellayout,int verbose){
  dvector Fy = newdvector(nsamples);
  double logL = multivariateregression(nvariables,nsamples,x,w,y,Fy,true,nullmodellayout,verbose);;
  freevector((void*)Fy);
  return logL;
}
