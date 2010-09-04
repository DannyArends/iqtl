/**
 * \file MYMATH/regression.cpp
 * \brief Code file, Implementation of: \ref LUdecomposition, \ref LUsolve, \ref LUinvert
 *
 * Copyright (c) 2010 Danny Arends
 * 
 **/
 
#include "regression.h"


double Lnormal(double residual, double variance){
  return exp(-pow(residual/sqrt(variance),2.0)/2.0 - log(sqrt(2.0*PI*variance)));
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

  backwardelimination(nvariables,nsamples,x,w,y);
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

double betacf(double a, double b, double x){
  double qap,qam,qab,em,tem,d,bz,bm=1.0,bp,bpp,az=1.0,am=1.0,ap,app,aold;
  int m;
  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  bz=1.0-qab*x/qap;
  for (m=1; m<=100; m++)
  {   em=(double)m;
      tem=em+em;
      d=em*(b-em)*x/((qam+tem)*(a+tem));
      ap=az+d*am;
      bp=bz+d*bm;
      d= -(a+em)*(qab+em)*x/((qap+tem)*(a+tem));
      app=ap+d*az;
      bpp=bp+d*bz;
      aold=az;
      am=ap/bpp;
      bm=bp/bpp;
      az=app/bpp;
      bz=1.0;
      if ( fabs((az-aold)/az)  < 3.0e-7) return az;
  }
  cout << "a or b too big or max number of iterations too small";
  exit(1); return 0.0;
}

/* functions gammln, betacf, betai necessary to calculate F(P,df1,df2) */
double gammln(double xx){
  double x,tmp,ser;
  static double cof[6]={76.18009173, -86.50532033, 24.01409822,
                     -1.231739516, 0.120858003e-2, -0.536382e-5};
  // if (xx<1) cout << "warning: full accuracy only for xx>1; xx= " << xx << endl;
  x=xx-1.0;
  tmp=x+5.5;
  tmp-= (x+0.5)*log(tmp);
  ser=1.0;
  for (int j=0; j<=5; j++) { x+=1.0; ser += cof[j]/x; }
  // delete[] cof;
  return -tmp+log(2.50662827465*ser);
}

double betai(double a, double b, double x){
  double bt;
  if (x<0.0 || x>1.0) { cout << "x not between 0 and 1"; exit(1); }
  if (x==0.0 || x==1.0) bt=0.0;
  else bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x<(a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
  else return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double inverseF(int df1, int df2, double alfa){
  double prob=0.0, minF=0.0, maxF=100.0, halfway=50.0, absdiff=1.0;
  int count=0;
  while ((absdiff>0.001)&&(count<100)){
    count++;
    halfway= (maxF+minF)/2.0;
    prob= betai(df2/2.0,df1/2.0,df2/(df2+df1*halfway));
    if (prob<alfa) maxF= halfway;
    else minF= halfway;
    absdiff= fabs(prob-alfa);
  }
  cout << "prob=" << prob << "; alfa=" << alfa << endl;
  return halfway;
}

void backwardelimination(uint nvariables,uint nsamples, dmatrix x, dvector w, dvector y){
  bool    finished   = false;
  uint    leastinterestingmodel;
  double  logLfull   = likelihoodbyem(nvariables,nsamples,x,w,y);
  bvector model      = newbvector(nvariables);
  double  dropneeded = 2*inverseF(2,nsamples-nvariables,0.005);
  cout << "Likelihood of the full model: " << logLfull << endl;
  while((!finished) && modelsize(nvariables,model) > 1){

    cout << "modelsize(model) = " << modelsize(nvariables,model) << "Drop " << dropneeded <<endl;
    dvector logL = newdvector(modelsize(nvariables,model));
    for(uint todrop=0;todrop<modelsize(nvariables,model);todrop++){
      bvector tempmodel = newbvector(nvariables);
      copybvector(nvariables,model,tempmodel);
      dropterm(nvariables,tempmodel,todrop);
      dmatrix designmatrix = createdesignmatrix(nvariables,nsamples,x,tempmodel);
      logL[todrop] = likelihoodbyem(modelsize(nvariables,tempmodel),nsamples,designmatrix,w,y);
      freematrix((void**)designmatrix,nsamples);
      freevector((void*)tempmodel);
    }

    leastinterestingmodel = lowestindex(modelsize(nvariables,model),logL);
    cout << "Least interesting model:" << leastinterestingmodel << " Difference to fullmodel:" << (logLfull - logL[leastinterestingmodel]) << endl;
    if(dropneeded > fabs(logLfull - logL[leastinterestingmodel])){
      dropterm(nvariables,model,leastinterestingmodel);
      logLfull = logL[leastinterestingmodel];
      cout << "Drop variable" << leastinterestingmodel << endl;
      cout << "Likelihood of the new full model: " << logLfull<< endl;
    }else{
      for(uint x=0;x<nvariables;x++){
        if(model[x]) cout << "Variable" << x << "In Model" << endl;
      }
      finished=true;
    }
  }

}

double likelihoodbyem(uint nvariables,uint nsamples, dmatrix x, dvector w, dvector y){
  uint   maxemcycles = 1000;
  uint   emcycle     = 0;
  double delta       = 1.0f;
  double logL        = 0.0f;
  double logLprev    = 0.0f;

  dvector Fy = newdvector(nsamples);
  //printdmatrix(x,nsamples,nvariables);
  while((emcycle<maxemcycles) && (delta>1.0e-5)){
    logL = multivariateregression(nvariables,nsamples,x,w,y,Fy);

    for(uint s=0;s<nsamples;s++){
      w[s] = (w[s]+Fy[s])/w[s];
    }

    delta= abs(logL-logLprev);
    logLprev=logL;
    emcycle++;
  }

  freevector((void*)Fy);

  cout << "[EM algorithm]\tFinished with "<< logL <<" after " << emcycle << "/" << maxemcycles << " cycles" << endl;
  return logL;
}

double multivariateregression(uint nvariables, uint nsamples, dmatrix x, dvector w, dvector y, dvector Fy){
  
  int d=0;
  double xtwj;
  dmatrix Xt   = newdmatrix(nvariables,nsamples);
  dmatrix XtWX = newdmatrix(nvariables, nvariables);
  dvector XtWY = newdvector(nvariables);

  ivector indx = newivector(nvariables);

  //cout << "calculating Xt" << endl;
  for(uint i=0; i<nsamples; i++){
    for(uint j=0; j<nvariables; j++){
      Xt[j][i] = x[i][j];
    }
  }

  //cout << "calculating XtWX and XtWY" << endl;
  for(uint i=0; i<nsamples; i++){
    for(uint j=0; j<nvariables; j++){
      xtwj     = Xt[j][i] * w[i];
      XtWY[j] += xtwj    * y[i];
      for(uint jj=0; jj<=j; jj++){
        XtWX[j][jj] += xtwj * Xt[jj][i];
      }
    }
  }
  
  LUdecomposition(XtWX, nvariables, indx, &d);
  LUsolve(XtWX, nvariables, indx, XtWY);

  //cout << "Estimated parameters:" << endl;
  //for (uint i=0; i < nvariables; i++){
  //  cout << "Parameter " << i << " = " << XtWY[i] << endl;
  //}

  dvector fit = newdvector(nsamples);
  dvector residual = newdvector(nsamples);
  dvector indL = newdvector(nsamples);
  
  double variance= 0.0;
  double logL=0.0;

  for (uint i=0; i<nsamples; i++){
    fit[i]= 0.0;
    for (uint j=0; j<nvariables; j++){
      fit[i]       += Xt[j][i] * XtWY[j];
      residual[i]   = y[i]-fit[i];
      variance     += w[i]*pow(residual[i],2.0);
    }
    Fy[i]     = Lnormal(residual[i],variance);
    indL[i]  += w[i]*Fy[i];
    logL     += log(indL[i]);
  }
  
  //cout << "Estimated response:" << endl;
  //printdvector(fit,nsamples);

  //cout << "Residuals:" << endl;
  //printdvector(residual,nsamples);

  //cout << "Estimated Fy:" << endl;
  //printdvector(Fy,nsamples);

  //cout << "Variance: " << variance << endl;
  //cout << "Loglikelihood: " << logL << endl;
  freematrix((void**)Xt,nvariables);
  freematrix((void**)XtWX, nvariables);
  freevector((void*)XtWY);
  freevector((void*)fit);
  freevector((void*)residual);
  freevector((void*)indL);
  return logL;
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
      cout << "Singular matrix" << endl;
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
      cout << "Singular matrix" << endl;
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
  freevector((void*)swap);
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
