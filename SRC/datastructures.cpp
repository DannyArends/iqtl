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
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <climits> 
#include <cfloat>

void *mycalloc(uint num, uint size) {
  void *buf;
  buf = calloc(num,size);
  if (buf) memset(buf,0,num*size);
  return buf;
}

bvector newbvector(uint dim) {
  bvector v = (bvector)calloc(dim,sizeof(bool));
  if (v==NULL) {
    cout << "Not enough memory for new vector of dimension" << (dim+1) << endl;
  }
  for(uint x=0;x<dim;x++){
    v[x]=true;
  }
  return v;
}

dvector newdvector(uint dim) {
  dvector v = (dvector)calloc(dim,sizeof(double));
  if (v==NULL) {
    cout << "Not enough memory for new vector of dimension" << (dim+1) << endl;
  }
  return v;
}

ivector newivector(uint dim) {
  ivector v = (ivector)calloc(dim,sizeof(int));
  if (v==NULL) {
    cout << "Not enough memory for new vector of dimension" << (dim+1) << endl;
  }
  return v;
}
cvector newcvector(uint dim) {
  cvector v = (cvector)calloc(dim,sizeof(char));
  if (v==NULL) {
    cout << "Not enough memory for new vector of dimension" << (dim+1) << endl;
  }
  return v;
}

void printdvector(dvector v, uint dim){
  for(uint i=0; i<dim; i++) {
    cout << v[i] << "\t";
  }
  cout << endl;
}
void printcvector(cvector v, uint dim){
  for(uint i=0; i<dim; i++) {
    cout << v[i] << "\t";
  }
  cout << endl;
}
void printivector(ivector v, uint dim){
  for(uint i=0; i<dim; i++) {
    cout << v[i] << "\t";
  }
  cout << endl;
}

dmatrix newdmatrix(uint rows, uint cols){
  dmatrix m = (dmatrix)calloc(rows,sizeof(dvector));
  if(m==NULL){
    cout << "Not enough memory for newmatrix" << endl;
  }
  for(uint i=0; i<rows; i++){
    m[i]= newdvector(cols);
  }
  return m;
}
cmatrix newcmatrix(uint rows, uint cols){
  cmatrix m = (cmatrix)calloc(rows,sizeof(cvector));
  if(m==NULL){
    cout << "Not enough memory for newmatrix" << endl;
  }
  for(uint i=0; i<rows; i++){
    m[i]= newcvector(cols);
  }
  return m;
}
imatrix newimatrix(uint rows, uint cols){
  imatrix m = (imatrix)calloc(rows,sizeof(ivector));
  if(m==NULL){
    cout << "Not enough memory for newmatrix" << endl;
  }
  for(uint i=0; i<rows; i++){
    m[i]= newivector(cols);
  }
  return m;
}

void printdmatrix(dmatrix m, uint rows, uint cols){
  for(uint r=0; r<rows; r++) {
    for(uint c=0; c<cols; c++) {
      cout << m[r][c] << "\t";
    }
    cout << endl;
  }
}
void printcmatrix(cmatrix m, uint rows, uint cols){
  for(uint r=0; r<rows; r++) {
    for(uint c=0; c<cols; c++) {
      cout << m[r][c] << "\t";
    }
    cout << endl;
  }
}
void printimatrix(imatrix m, uint rows, uint cols){
  for(uint r=0; r<rows; r++) {
    for(uint c=0; c<cols; c++) {
      cout << m[r][c] << "\t";
    }
    cout << endl;
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
