/**
 * \file expectationmaximization.cpp
 * \brief Code file containing functions
 *
 * last modified Sep, 2010
 * first written Apr, 2010
 * Copyright (c) 2010 Danny Arends
 * 
 **/
 
#include "expectationmaximization.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <limits> 
#include <cfloat>
#include <R.h>
#include <Rmath.h>

void nullmodellikelihood_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* nullmodellayout,int* verbose,double* out){
  dmatrix transformedx;
  dvectortodmatrix((*nsamples),(*nvariables),x,&transformedx);
 (*out) = 2 * nullmodel((*nvariables), (*nsamples), transformedx, w, y, nullmodellayout, (*verbose));
 if((*verbose)) Rprintf("null likelihood: %f\n",(*out));

}

void modellikelihoodbyem_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* verbose,double* out){
  dmatrix transformedx;
  dvectortodmatrix((*nsamples),(*nvariables),x,&transformedx);
 (*out) = 2 * likelihoodbyem((*nvariables), (*nsamples), transformedx, w, y, (*verbose));
 if((*verbose)) Rprintf("model likelihood: %f\n",(*out));

}

void lodscorebyem_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* nullmodellayout,int* verbose,double* out){
  dmatrix transformedx;
  dvectortodmatrix((*nsamples),(*nvariables),x,&transformedx);
  (*out) = (2*likelihoodbyem((*nvariables), (*nsamples), transformedx, w, y, (*verbose)) - 2 * nullmodel((*nvariables), (*nsamples), transformedx, w, y, nullmodellayout, (*verbose))) / 4.60517;
  if((*verbose)) Rprintf("lodscore: %f\n",(*out));

}

double likelihoodbyem(uint nvariables,uint nsamples, dmatrix x, dvector w, dvector y,int verbose){
  uint   maxemcycles = 1000;
  uint   emcycle     = 0;
  double delta       = 1.0f;
  double logL        = 0.0f;
  double logLprev    = 0.0f;
  
  if(verbose){
    Rprintf("DesignMatrix:\n");
    printdmatrix(x,nsamples,nvariables);
  }
  dvector Fy = newdvector(nsamples);
  ivector nullmodellayout = newivector(nvariables);
  if(verbose) Rprintf("Starting EM\n");
  while((emcycle<maxemcycles) && (delta > 1.0e-9)){
    logL = multivariateregression(nvariables,nsamples,x,w,y,Fy,false,nullmodellayout,verbose);
    for(uint s=0;s<nsamples;s++){
      if(w[s] != 0) w[s] = (w[s] + Fy[s])/w[s];
    }
    delta = fabs(logL-logLprev);
    logLprev=logL;
    emcycle++;
  }
  freevector((void*)Fy);
  
  //Rprintf("EM %d/%d cycles\n",emcycle,maxemcycles);
  return (logL);
}
