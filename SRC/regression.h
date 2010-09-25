/**
 * \file regression.h
 * \brief Header file, for regression.cpp
 *
 * last modified Sep, 2010
 * first written Apr, 2010 
 * Copyright (c) 2010 Danny Arends
 * 
 **/

 #ifndef REGRESSION_H_
  #define REGRESSION_H_
extern "C"
{  
  #include "datastructures.h"

  double multivariateregression(uint nvariables, uint nsamples, dmatrix x, dvector w, dvector y, dvector Fy,bool nullmodel,ivector nullmodellayout,int verbose);
  double nullmodel(uint nvariables, uint nsamples, dmatrix x, dvector w, dvector y,ivector nullmodellayout,int verbose);
}
#endif
