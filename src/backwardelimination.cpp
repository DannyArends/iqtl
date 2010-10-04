/**
 * \file backwardelimination.cpp
 * \brief Code file containing functions
 *
 * last modified Sep, 2010
 * first written Apr, 2010
 * Copyright (c) 2010 Danny Arends
 * 
 **/

#include "backwardelimination.h" 
#include "regression.h"
#include "expectationmaximization.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <limits> 
#include <cfloat>
#include <R.h>
#include <Rmath.h>

void testregression(){
  uint   nvariables  = 10;
  uint   nsamples    = 300;
  dmatrix x = newdmatrix(nsamples,nvariables);
  dvector w = newdvector(nsamples);
  dvector y = newdvector(nsamples);
  for(uint s=0;s<nsamples;s++){
    for(uint v=0;v<nvariables;v++){
      x[s][v] = float(rand()%2)-1.0f/2.0f;
    }
    w[s]= 1.0f;
    y[s]= s + float(rand()%10);
  }

  backwardelimination(nvariables,nsamples,x,w,y,1);
}

uint modelsize(uint nvariables, bvector model){
  uint s=0;
  for(uint v=0;v<nvariables;v++){
    if(model[v]) s++;
  }
  return s;
}

void dropterm(uint nvariables, bvector model,uint which){
  for(uint v=0;v<nvariables;v++){
    if(model[v]){
      if(which==0){
        model[v]=false;
        return;
      }else{
        which--;
      }
    }
  }
}

uint lowestindex(uint dim, dvector values){
  double min = values[0];
  uint index = 0;
  for(uint v=0;v<dim;v++){
    if(values[v] < min){
      min=values[v];
      index=v;
    }
  }
  return index;
}

uint highestindex(uint dim, dvector values){
  double max = values[0];
  uint index = 0;
  for(uint v=0;v<dim;v++){
    if(values[v] > max){
      max=values[v];
      index=v;
    }
  }
  return index;
}

void copybvector(uint dim,bvector origin,bvector target){
  for(uint v=0;v<dim;v++){
    target[v]=origin[v];
  }
}

dmatrix createdesignmatrix(uint nvariables,uint nsamples, dmatrix x, bvector model){
  uint    dimension = modelsize(nvariables,model);
  dmatrix Xn = newdmatrix(nsamples,dimension);
  uint    newcol=0;

  for(uint v=0;v<nvariables;v++){
    if(model[v]){
      for(uint s=0;s<nsamples;s++){
        Xn[s][newcol] = x[s][v];
      }
      newcol++;
    }
  }
  return Xn;
}

void inverseF_R(int* df1,int* df2, double* alfa, double* out){
  (*out) = inverseF((*df1), (*df1), (*alfa));
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

void backwardelimination_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* verbose,double* out){
  dmatrix transformedx;
  dvectortodmatrix((*nsamples),(*nvariables),x,&transformedx);
  backwardelimination((*nvariables), (*nsamples), transformedx, w, y, (*verbose));
}

void backwardelimination(uint nvariables,uint nsamples, dmatrix x, dvector w, dvector y,int verbose){
  bool    finished   = false;
  uint    leastinterestingmodel;
  double  logLfull   = likelihoodbyem(nvariables,nsamples,x,w,y,0);
  bvector model      = newbvector(nvariables);
  double  dropneeded = 2*inverseF(2,nsamples-nvariables,0.005);
  if(verbose) Rprintf("Likelihood of the full model: %f\n",logLfull);
  while((!finished) && modelsize(nvariables,model) > 1){

    if(verbose) Rprintf("modelsize(model) = %d, Drop %f\n", modelsize(nvariables,model),dropneeded);
    dvector logL = newdvector(modelsize(nvariables,model));
    for(uint todrop=0;todrop<modelsize(nvariables,model);todrop++){
      bvector tempmodel = newbvector(nvariables);
      copybvector(nvariables,model,tempmodel);
      dropterm(nvariables,tempmodel,todrop);
      dmatrix designmatrix = createdesignmatrix(nvariables,nsamples,x,tempmodel);
      logL[todrop] = likelihoodbyem(modelsize(nvariables,tempmodel),nsamples,designmatrix,w,y,0);
      freematrix((void**)designmatrix,nsamples);
      freevector((void*)tempmodel);
    }

    leastinterestingmodel = highestindex(modelsize(nvariables,model),logL);
    if(verbose) Rprintf("Least interesting model: %d\n", leastinterestingmodel);
    if(verbose) Rprintf("Difference to fullmodel: %f\n", (logLfull - logL[leastinterestingmodel]));
    if(dropneeded > fabs(logLfull - logL[leastinterestingmodel])){
      dropterm(nvariables,model,leastinterestingmodel);
      logLfull = logL[leastinterestingmodel];
      if(verbose) Rprintf("Drop variable %d\n", leastinterestingmodel);
      if(verbose) Rprintf("Likelihood of the new full model: %f\n",logLfull);
    }else{
      Rprintf("\n\nWe have a model\n",x);
      for(uint x=0;x<nvariables;x++){
        if(model[x]) Rprintf("Variable %d in Model\n",x);
      }
      finished=true;
    }
  }
}

