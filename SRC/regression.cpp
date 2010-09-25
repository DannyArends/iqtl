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
      if(nullmodellayout[(i-1)] == 1){
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

double nullmodel(uint nvariables, uint nsamples, dmatrix x, dvector w, dvector y,ivector nullmodellayout,int verbose){
  dvector Fy = newdvector(nsamples);
  double logL = multivariateregression(nvariables,nsamples,x,w,y,Fy,true,nullmodellayout,verbose);;
  freevector((void*)Fy);
  return logL;
}
