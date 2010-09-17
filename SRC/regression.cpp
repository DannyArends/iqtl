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
  Rprintf("prob = %f, alfa= %f\n",prob,alfa);
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
      if(verbose) Rprintf("Likelihood of the new full model: %f",logLfull);
    }else{
      for(uint x=0;x<nvariables;x++){
        if(model[x]) Rprintf("Variable %d in Model",x);
      }
      finished=true;
    }
  }
}

void modellikelihoodbyem_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* verbose,double* out){
  dmatrix transformedx;
  dvectortodmatrix((*nsamples),(*nvariables),x,&transformedx);
 (*out) = likelihoodbyem((*nvariables), (*nsamples), transformedx, w, y, (*verbose));
}

void nulllikelihoodbyem_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* verbose,double* out){
  dmatrix transformedx;
  dvectortodmatrix((*nsamples),(*nvariables),x,&transformedx);
 (*out) = nullmodel((*nvariables), (*nsamples), transformedx, w, y, (*verbose));
}

void likelihoodbyem_R(int* nvariables,int* nsamples, double* x, double* w, double* y,int* verbose,double* out){
  dmatrix transformedx;
  dvectortodmatrix((*nsamples),(*nvariables),x,&transformedx);
 (*out) = likelihoodbyem((*nvariables), (*nsamples), transformedx, w, y, (*verbose));
 (*out) -= nullmodel((*nvariables), (*nsamples), transformedx, w, y, (*verbose));
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
  
  if(verbose) Rprintf("Starting EM\n");
  while((emcycle<maxemcycles) && (delta > 1.0e-9)){
    logL = multivariateregression(nvariables,nsamples,x,w,y,Fy,verbose);
    for(uint s=0;s<nsamples;s++){
      if(w[s] != 0) w[s] = (w[s] + Fy[s])/w[s];
    }
    delta = fabs(logL-logLprev);
    logLprev=logL;
    emcycle++;
  }
  freevector((void*)Fy);
  
  Rprintf("Finished with %f after %d/%d cycles\n",logL,emcycle,maxemcycles);
  return (logL);
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

double multivariateregression(uint nvariables, uint nsamples, dmatrix x, dvector w, dvector y, dvector Fy,int verbose){
  dmatrix Xt   = translatematrix(nvariables,nsamples,x,verbose);
  if(verbose){
    Rprintf("T(DesignMatrix):\n");
    printdmatrix(Xt,nvariables,nsamples);
  }
  dvector XtWY = calculateparameters(nvariables,nsamples,Xt,w,y,verbose);

  if(verbose){
    Rprintf("Estimated parameters:\n");
    printdvector(XtWY,nvariables);
  }

  double variance= 0.0;
  double logLQTL=0.0;
  dvector fit       = newdvector(nsamples);
  dvector residual  = newdvector(nsamples);
  ldvector indL     = newldvector(nsamples);
  
  for (uint i=0; i<nsamples; i++){
    fit[i]= 0.0;
    for (uint j=0; j<nvariables; j++){
      fit[i]       += Xt[j][i] * XtWY[j];
      residual[i]   = y[i]-fit[i];
      variance     += w[i]*pow(residual[i],2.0);
    }
    Fy[i]     = Lnormal(residual[i],variance);
    indL[i]  += w[i]*Fy[i];
    logLQTL  += log10(indL[i]);
  }
  
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
  freevector((void*)indL);

  return logLQTL;
}

double nullmodel(uint nvariables, uint nsamples, dmatrix x, dvector w, dvector y,int verbose){
  dmatrix Xt   = translatematrix(nvariables,nsamples,x,verbose);
  dvector XtWY = calculateparameters(nvariables,nsamples,Xt,w,y,verbose);
  dvector Fy   = newdvector(nsamples);
  
  if(verbose){
    Rprintf("Estimated parameters:\n");
    printdvector(XtWY,nvariables);
  }

  double variance= 0.0;
  double logLNULL=0.0;
  dvector fit = newdvector(nsamples);
  dvector residual = newdvector(nsamples);
  ldvector indL = newldvector(nsamples);

  for (uint i=1; i < nvariables; i++){
    XtWY[i] = 0.0;
  }

  for (uint i=0; i<nsamples; i++){
    fit[i]= 0.0;
    for (uint j=0; j<nvariables; j++){
      fit[i]       += Xt[j][i] * XtWY[j];
      residual[i]   = y[i]-fit[i];
      variance     += w[i]*pow(residual[i],2.0);
    }
    Fy[i]     = Lnormal(residual[i],variance);
    indL[i]  += w[i]*Fy[i];
    logLNULL += log10(indL[i]);
  }
  
  if(verbose){
    Rprintf("Estimated response:\n");
    printdvector(fit,nsamples);

    Rprintf("Residuals:\n");
    printdvector(residual,nsamples);

    Rprintf("Estimated Fy:\n");
    printdvector(Fy,nsamples);

    Rprintf("Variance: %f\n",variance);
    Rprintf("Loglikelihood NULL model: %f\n",logLNULL);
  }
  
  freematrix((void**)Xt,nvariables);
  freevector((void*)XtWY);
  freevector((void*)fit);
  freevector((void*)Fy);
  freevector((void*)residual);
  freevector((void*)indL);
  
  return logLNULL;
}

bool LUdecomposition(dmatrix m, int dim, ivector ndx, int *d) {
  int r, c, rowmax, i;
  double max, temp, sum;
  dvector swap = newdvector(dim);
  dvector scale = newdvector(dim);
  *d=1;
  for (r=0; r<dim; r++) {
    for (max=0.0, c=0; c<dim; c++){
      if ((temp=fabs(m[r][c])) > max){
        max=temp;
      }
    }
    if (max==0.0){
      Rprintf("Singular matrix\n");
      return false;
    }
    scale[r]=1.0/max;
  }
  for (c=0; c<dim; c++) {
    for (r=0; r<c; r++) {
      for (sum=m[r][c], i=0; i<r; i++) sum-= m[r][i]*m[i][c];
      m[r][c]=sum;
    }
    for (max=0.0, rowmax=c, r=c; r<dim; r++) {
      for (sum=m[r][c], i=0; i<c; i++) sum-= m[r][i]*m[i][c];
      m[r][c]=sum;
      if ((temp=scale[r]*fabs(sum)) > max) {
        max=temp;
        rowmax=r;
      }
    }
    if (max==0.0){
      Rprintf("Singular matrix\n");
      return false;
    }
    if (rowmax!=c) {
      swap=m[rowmax];
      m[rowmax]=m[c];
      m[c]=swap;
      scale[rowmax]=scale[c];
      (*d)= -(*d);
    }
    ndx[c]=rowmax;
    temp=1.0/m[c][c];
    for(r=c+1; r<dim; r++){
      m[r][c]*=temp;
    }
  }
  freevector((void*)scale);
  return true;
}

void LUsolve(dmatrix lu, int dim, ivector ndx, dvector b) {
  int r, c;
  double sum;
  for (r=0; r<dim; r++) {
    sum=b[ndx[r]];
    b[ndx[r]]=b[r];
    for (c=0; c<r; c++) sum-= lu[r][c]*b[c];
    b[r]=sum;
  }
  for (r=dim-1; r>-1; r--) {
    sum=b[r];
    for (c=r+1; c<dim; c++) sum-= lu[r][c]*b[c];
    b[r]=sum/lu[r][r];
  }
}

void LUinvert(dmatrix lu, dmatrix inv, int dim, int *ndx){
  int r,c;
  dvector b = newdvector(dim);
  for (c=0; c<dim; c++){
     b[c]=1.0;
     LUsolve(lu,dim,ndx,b);
     for (r=0; r<dim; r++) inv[r][c]= b[r];
  }
  freevector((void*)b);
} 
