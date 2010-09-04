/**********************************************************************
 *
 * datastructures.cpp
 *
 * Copyright (c) 2010 Danny Arends
 *
 **********************************************************************/
#include "datastructures.h"

void *mycalloc(uint num, uint size) {
  void *buf;
  buf = calloc(num,size);
  if (buf) memset(buf,0,num*size);
  return buf;
}

bvector newbvector(uint dim) {
  bvector v = (bvector)calloc(dim,sizeof(bool));
  if (v==NULL) {
    Rprintf("Not enough memory for new vector of dimension %i\n",dim+1);
  }
  for(uint x=0;x<dim;x++){
    v[x]=true;
  }
  return v;
}

dvector newdvector(uint dim) {
  dvector v = (dvector)calloc(dim,sizeof(double));
  if (v==NULL) {
    Rprintf("Not enough memory for new vector of dimension %i\n",dim+1);
  }
  return v;
}

ivector newivector(uint dim) {
  ivector v = (ivector)calloc(dim,sizeof(int));
  if (v==NULL) {
    Rprintf("Not enough memory for new vector of dimension %i\n",dim+1);
  }
  return v;
}
cvector newcvector(uint dim) {
  cvector v = (cvector)calloc(dim,sizeof(char));
  if (v==NULL) {
    Rprintf("Not enough memory for new vector of dimension %i\n",dim+1);
  }
  return v;
}

void printdvector(dvector v, uint dim){
  for(uint i=0; i<dim; i++) {
    Rprintf("%f\t",v[i]);
  }
  Rprintf("\n");
}
void printcvector(cvector v, uint dim){
  for(uint i=0; i<dim; i++) {
    Rprintf("%c\t",v[i]);
  }
  Rprintf("\n");
}
void printivector(ivector v, uint dim){
  for(uint i=0; i<dim; i++) {
    Rprintf("%i\t",v[i]);
  }
  Rprintf("\n");
}

dmatrix newdmatrix(uint rows, uint cols){
  dmatrix m = (dmatrix)calloc(rows,sizeof(dvector));
  if(m==NULL){
    Rprintf("Not enough memory for new matrix\n");
  }
  for(uint i=0; i<rows; i++){
    m[i]= newdvector(cols);
  }
  return m;
}
cmatrix newcmatrix(uint rows, uint cols){
  cmatrix m = (cmatrix)calloc(rows,sizeof(cvector));
  if(m==NULL){
    Rprintf("Not enough memory for new matrix\n");
  }
  for(uint i=0; i<rows; i++){
    m[i]= newcvector(cols);
  }
  return m;
}
imatrix newimatrix(uint rows, uint cols){
  imatrix m = (imatrix)calloc(rows,sizeof(ivector));
  if(m==NULL){
    Rprintf("Not enough memory for new matrix\n");
  }
  for(uint i=0; i<rows; i++){
    m[i]= newivector(cols);
  }
  return m;
}

void printdmatrix(dmatrix m, uint rows, uint cols){
  for(uint r=0; r<rows; r++) {
    for(uint c=0; c<cols; c++) {
      Rprintf("%f\t",m[r][c]);
    }
    Rprintf("\n");
  }
}
void printcmatrix(cmatrix m, uint rows, uint cols){
  for(uint r=0; r<rows; r++) {
    for(uint c=0; c<cols; c++) {
      Rprintf("%c\t",m[r][c]);
    }
    Rprintf("\n");
  }
}
void printimatrix(imatrix m, uint rows, uint cols){
  for(uint r=0; r<rows; r++) {
    for(uint c=0; c<cols; c++) {
      Rprintf("%i\t",m[r][c]);
    }
    Rprintf("\n");
  }
}

void freevector(void *v) {
  free(v);
}
void freematrix(void **m,uint rows) {
  for (uint i=0; i < rows; i++) {
    free(m[i]);
  }
  free(m);
}