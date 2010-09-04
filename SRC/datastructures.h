/**********************************************************************
 *
 * mqmdatatypes.h
 *
 * Copyright (c) 2010 Danny Arends
 *
 **********************************************************************/
#ifndef _DATASTRUCTS_H
  #define _DATASTRUCTS_H
extern "C"
{  
  #include <R.h>
  typedef unsigned int  uint;
  typedef double*       dvector;
  typedef dvector*      dmatrix;

  dvector newdvector(uint dim);
  dmatrix newdmatrix(uint rows, uint cols);
  void    printdmatrix(dmatrix m, uint rows, uint cols);
  void    freematrix(void **m, uint rows);
  void    freedmatrix(dmatrix m, uint rows);
  void    *mycalloc(uint num, uint size);
}
#endif
