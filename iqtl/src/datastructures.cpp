/**
 * \file datastructures.cpp
 * \brief Code file containing functions to manipulate datastructures
 *
 * last modified Sep, 2010
 * first written Apr, 2010
 * copyright (c) 2010 Danny Arends
 *
 * C functions for the iqtl package
 * Contains:  newbvector,newdvector,newcvector,newivector
 *            printdvector,printcvector,printivector
 *            newdmatrix,newcmatrix,newimatrix
 *            printdmatrix,printcmatrix,printimatrix
 *            freematrix,freevector,mycalloc
 *
 **/
 
#include "datastructures.h"
#include <assert.h>
#include <cstdlib>
#include <fstream>
#include <climits> 
#include <cfloat>

void *mycalloc(uint num, uint size) {
  void *buf;
  buf = R_chk_calloc(num,size);
  if (buf) memset(buf,0,num*size);
  return buf;
}

bvector newbvector(uint dim) {
  bvector v = (bvector)mycalloc(dim,sizeof(bool));
  if(v==NULL){ Rprintf("Not enough memory for new vector of dimension %d\n", (dim+1)); }
  for(size_t x = 0; x < dim; x++){ v[x]=true; }
  return v;
}

ldvector newldvector(uint dim) {
  ldvector v = (ldvector)mycalloc(dim,sizeof(long double));
  if(v==NULL){ Rprintf("Not enough memory for new vector of dimension %d\n", (dim+1)); }
  return v;
}

dvector newdvector(uint dim) {
  dvector v = (dvector)mycalloc(dim,sizeof(double));
  if(v==NULL){ Rprintf("Not enough memory for new vector of dimension %d\n", (dim+1)); }
  return v;
}

ivector newivector(uint dim) {
  ivector v = (ivector)mycalloc(dim,sizeof(int));
  if(v==NULL){ Rprintf("Not enough memory for new vector of dimension %d\n", (dim+1)); }
  return v;
}

cvector newcvector(uint dim) {
  cvector v = (cvector)mycalloc(dim,sizeof(char));
  if(v==NULL){ Rprintf("Not enough memory for new vector of dimension %d\n", (dim+1)); }
  return v;
}

void printdvector(dvector v, uint dim){
  for(size_t i=0; i<dim; i++){ Rprintf("%f\t",v[i]); }
  Rprintf("\n");
}

void printldvector(dvector v, uint dim){
  for(size_t i=0; i<dim; i++){ Rprintf("%f\t",v[i]); }
  Rprintf("\n");
}

void printcvector(cvector v, uint dim){
  for(size_t i=0; i<dim; i++){ Rprintf("%c\t",v[i]); }
  Rprintf("\n");
}

void printivector(ivector v, uint dim){
  for(size_t i=0; i<dim; i++){ Rprintf("%d\t",v[i]); }
  Rprintf("\n");
}

dmatrix newdmatrix(uint rows, uint cols){
  dmatrix m = (dmatrix)mycalloc(rows,sizeof(dvector));
  if(m==NULL){ Rprintf("Not enough memory for a new matrix\n"); }
  for(size_t i = 0; i < rows; i++){ m[i]= newdvector(cols); }
  return m;
}

cmatrix newcmatrix(uint rows, uint cols){
  cmatrix m = (cmatrix)mycalloc(rows,sizeof(cvector));
  if(m==NULL){ Rprintf("Not enough memory for new matrix\n"); }
  for(size_t i = 0; i < rows; i++){ m[i]= newcvector(cols); }
  return m;
}

imatrix newimatrix(uint rows, uint cols){
  imatrix m = (imatrix)mycalloc(rows, sizeof(ivector));
  if(m==NULL){ Rprintf("Not enough memory for new matrix\n"); }
  for(size_t i = 0; i < rows; i++){ m[i]= newivector(cols); }
  return m;
}

void printdmatrix(dmatrix m, uint rows, uint cols){
  for(size_t r = 0; r<rows; r++) {
    for(size_t c = 0; c<cols; c++) {
      Rprintf("%f\t",m[r][c]);
    }
    Rprintf("\n");
  }
}
void printcmatrix(cmatrix m, uint rows, uint cols){
  for(size_t r=0; r < rows; r++){
    for(size_t c=0; c < cols; c++){ Rprintf("%c\t", m[r][c]); }
    Rprintf("\n");
  }
}

void printimatrix(imatrix m, uint rows, uint cols){
  for(size_t r=0; r < rows; r++){
    for(size_t c=0; c < cols; c++){ Rprintf("%d\t", m[r][c]); }
    Rprintf("\n");
  }
}

void dvectortodmatrix(int nrow, int ncol, dvector in, dmatrix* out){
  *out = (double**)R_alloc(nrow, sizeof(double*));
  for(size_t r = 0; r < nrow; r++){
    (*out)[r] = (double*)R_alloc(ncol, sizeof(double));
    for(size_t c = 0; c < ncol; c++){
      (*out)[r][c] = in[(c*nrow)+r];
    }
  }
}

void ivectortoimatrix(int nrow, int ncol, ivector in, imatrix* out){
  *out = (int**)R_alloc(nrow, sizeof(int*));
  for(size_t r = 0; r < nrow; r++){
    (*out)[r] = (int*)R_alloc(ncol, sizeof(int));
    for(size_t c = 0; c < ncol; c++){
      (*out)[r][c] = in[(c*nrow)+r];
    }
  }
}

void freevector(void *v){ if(v != NULL) Free(v); }

void freematrix(void **m,uint rows) {
  for(size_t i = 0; i < rows; i++){ freevector(m[i]); }
  if(m != NULL) Free(m);
}
