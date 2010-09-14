/**********************************************************************
 * 
 * xgapbin_r.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Apr, 2009
 * first written Apr, 2009
 *
 * C functions for the iqtl package
 * Contains: R_load_XGAPheader
 *           R_load_XGAPnames
 *           R_load_XGAPdouble
 *           R_load_XGAPstring
 *
 **********************************************************************/
 
 #ifndef _XGAPBIN_R_H
  #define _XGAPBIN_R_H
  extern "C"
  {

  #include "datastructures.h"
  #include "xgapparser.h"

  void R_load_XGAPheader(char **filename,int* num_rows,int* num_cols,int* isNumeric,int* pos);
  void R_load_XGAPnames (char **filename,int* num_rows,int* num_cols,int* pos, char** rownames, char** colnames);
  void R_load_XGAPdouble(char **filename,int* num_rows,int* num_cols,int* pos, double* data);
  void R_load_XGAPstring(char **filename,int* num_rows,int* num_cols,int* pos, char** data);

  }
#endif
