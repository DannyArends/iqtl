/**
 * \file LUdecomposition.h
 * \brief Header file, for LUdecomposition.cpp
 *
 * last modified Sep, 2010
 * first written Apr, 2010 
 * Copyright (c) 2010 Danny Arends
 * 
 **/

 #ifndef LUDECOMPOSITION_H_
  #define LUDECOMPOSITION_H_
extern "C"
{  
  #include "datastructures.h"

  bool LUdecomposition(dmatrix m, int dim, ivector ndx, int *d);
  void LUsolve(dmatrix lu, int dim, ivector ndx, dvector b);
  void LUinvert(dmatrix lu, dmatrix inv, int dim, int *ndx);
}
#endif
