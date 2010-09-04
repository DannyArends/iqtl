/**********************************************************************
 * 
 * R_interface.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Apr, 2009
 * first written Apr, 2009
 *
 * C functions for the Rexamples package
 * Contains: R_add_in_C
 *
 **********************************************************************/
extern "C"
{
#include <R.h>
#include <Rinternals.h>
#include "datastructures.h"
#include "parser.h"

void R_load_XGAPheader(char **filename,int* num_rows,int* num_cols,int* isNumeric,int* pos);
void R_load_XGAPnames (char **filename,int* num_rows,int* num_cols,int* pos, char** rownames, char** colnames);
void R_load_XGAPdouble(char **filename,int* num_rows,int* num_cols,int* pos, double* data);
void R_load_XGAPstring(char **filename,int* num_rows,int* num_cols,int* pos, char** data);

}
