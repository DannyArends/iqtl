/**********************************************************************
 *
 * mqmdatatypes.h
 *
 * Copyright (c) 2010 Danny Arends
 *
 **********************************************************************/
#ifndef _DATASTRUCTS_H
  #define _DATASTRUCTS_H
  #include "typedefs.h"
 extern "C"
{
bvector newbvector(uint dim);
dvector newdvector(uint dim);
cvector newcvector(uint dim);
ivector newivector(uint dim);

void printdvector(dvector v, uint dim);
void printcvector(cvector v, uint dim);
void printivector(ivector v, uint dim);

dmatrix newdmatrix(uint rows, uint cols);
cmatrix newcmatrix(uint rows, uint cols);
imatrix newimatrix(uint rows, uint cols);

void printdmatrix(dmatrix m, uint rows, uint cols);
void printcmatrix(cmatrix m, uint rows, uint cols);
void printimatrix(imatrix m, uint rows, uint cols);

void freematrix(void **m, uint rows);
void freevector(void *v);
void *mycalloc(uint num, uint size);
}
#endif
