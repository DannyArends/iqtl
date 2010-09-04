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

dvector newdvector(uint dim){
  dvector v;
  v = (dvector)mycalloc(dim, sizeof(double));
  if (v==NULL) {
    Rprintf("Not enough memory for new vector of dimension %d\n",(dim+1));
  }
  return v;
}

dmatrix newdmatrix(uint rows, uint cols) {
  dmatrix m;
  m = (dmatrix)mycalloc(rows, sizeof(dvector));
  if (m==NULL) {
    Rprintf("Not enough memory for new vector of dimension %d %d\n",(rows+1),(cols+1));
  }
  for (uint i=0; i<rows; i++) {
    m[i]= newdvector(cols);
  }
  return m;
}

void printdmatrix(dmatrix m, uint rows, uint cols) {
  for (uint r=0; r<rows; r++) {
    for (uint c=0; c<cols; c++) {
      //Rprintf("%f\t",m[r][c]);
    }
    //Rprintf("\n");
  }
}

void freevector(void *v) {
  free(v);
}
void freematrix(void **m, uint rows) {
  for (uint i=0; i<rows; i++) {
    free(m[i]);
  }
  free(m);
}

void freedmatrix(dmatrix m, uint rows) {
  freematrix((void**)m,rows);
}
