/**
 * \file xgapwriter.cpp
 * \brief Functions to parse a XGAP binary file
 *
 * last modified Jun, 2011
 * first written Jun, 2011
 * copyright (c) 2009-2010 Danny Arends
 *
 * C functions for the iqtl package
 * Contains: xgapwriter
 *
 **/
 
#include "xgapwriter.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits> 
#include <cfloat>

enum DataType{ EMPTY = 0, INTMATRIX = 1, DOUBLEMATRIX = 2, FIXEDCHARMATRIX = 3, VARCHARMATRIX = 4, MATRIXREFERENCE = 5};
uint     footprint[2] = { 0, 5 };
uint     sprint[2] = { 1, 1 };
uint     eprint[2] = { 0, 0 };
int      version[3] = { 0, 0, 1 };

struct matrixHeader{
  int namelength;
  const char* name;
  
  DataType type;
  int columns;
  int rows;
};


struct intMatrix{
  int** data;
};

struct doubleMatrix{
  double** data;
};

struct fixedcharMatrix{
  char** data;
  int length;
};

struct varcharMatrix{
  char*** data;
  int lengths[];
};

struct refMatrix{
  int** colpointers;
};

void writeMatrixHeader(ofstream* myfile, matrixHeader h){
  (*myfile).write((char*) &h.namelength,  sizeof(int));
  (*myfile).write(         h.name,        sizeof(char) * h.namelength);
  (*myfile).write((char*) &h.type,        sizeof(DataType));
  (*myfile).write((char*) &h.columns,     sizeof(int));
  (*myfile).write((char*) &h.rows,        sizeof(int));
}

void writeintMatrixData(ofstream* myfile, matrixHeader h, intMatrix matrix){
  for(int r;r<h.rows;r++){
    (*myfile).write((char*) &matrix.data[r], h.columns* sizeof(int));
  }
}

void writedoubleMatrixData(ofstream* myfile, matrixHeader h, doubleMatrix matrix){
  for(int r;r<h.rows;r++){
    (*myfile).write((char*) &matrix.data[r], h.columns* sizeof(double));
  }
}

void writefixedCharMatrixData(ofstream* myfile, matrixHeader h, fixedcharMatrix matrix){
  (*myfile).write((char*) &matrix.length, sizeof(int));
  for(int r;r<h.rows;r++){
    for(int c;c<h.columns;c++){  
      (*myfile).write((char*) &matrix.data[r][c], matrix.length * sizeof(char));
    }
  }
}

void writevarCharMatrixData(ofstream* myfile, matrixHeader h, varcharMatrix matrix){
  (*myfile).write((char*) &matrix.lengths, h.rows*h.columns* sizeof(int));
  for(int r;r<h.rows;r++){
    for(int c;c<h.columns;c++){  
      (*myfile).write((char*) &matrix.data[r][c], matrix.lengths[r+c*h.rows]*sizeof(char));
    }
  }
}


void writeReferenceMatrixData(ofstream* myfile, matrixHeader h, void* data){

}


void reorg_data(int rows, int cols, int type,int* lengths, void* in, int*** integers, double*** doubles, char*** text) {
//reorganisation of integers into a matrix
  int i;
  if(type == INTMATRIX){ 
    *integers = (int **)R_alloc(cols, sizeof(int*));
    (*integers)[0] = (int*)in;
    for (i=1; i< cols; i++)
      (*integers)[i] = (*integers)[i-1] + rows;
  }
  if(type == DOUBLEMATRIX){ 
    *doubles = (double **)R_alloc(cols, sizeof(double*));
    (*doubles)[0] = (double*)in;
    for (i=1; i< cols; i++)
      (*doubles)[i] = (*doubles)[i-1] + rows;
  }
  if(type == FIXEDCHARMATRIX){ 
    *text = (char **)R_alloc(cols, (*lengths) * sizeof(char*));
    (*text)[0] = (char*)in;
    for (i=1; i< cols; i++)
      (*text)[i] = (*text)[i-1] + rows;
  }
}

void r_writebinaryfile(char** filename, int* namelength, int* nrows, int* ncols, int* type,int** lengths, void** data){
  Rprintf("INFO: Starting C-part of the writebinaryfile\n");
  int** intdata;
  double** doubledata;
  char** chardata;
  reorg_data((*nrows),(*ncols),(*type),(*lengths),(*data),&intdata,&doubledata,&chardata);
  for(int c=0;c<5;c++){
    for(int r=0;r<5;r++){
      if((*type) == INTMATRIX) Rprintf("integer[%d,%d] = %d\n",r,c,intdata[r][c]);
      if((*type) == DOUBLEMATRIX) Rprintf("double[%d,%d] = %f\n",r,c,doubledata[r][c]);
      if((*type) == FIXEDCHARMATRIX) Rprintf("char[%d,%d] = %s\n",r,c,chardata[r][c]);
    }
  }
  ofstream myfile;
  myfile.open(filename[0], std::ios::out | std::ios::binary |std::ios::ate );
  if(myfile.is_open()){
    myfile.write((char*) &footprint,   sizeof(uint)*2);
    myfile.write((char*) &version,     sizeof(uint)*3);
    myfile.write((char*)  namelength,  sizeof(int));
    myfile.write(         filename[0], sizeof(char) * (*namelength));
    myfile.write((char*) &footprint,   sizeof(uint)*2);
    myfile.close();
  }else{
    Rprintf("Error opening file\n");
  }
  Rprintf("INFO: Done in c\n");
}

bool writeFile(const char* filename, int namelength, int nmatrices, matrixHeader* headers, void** data){
  //cout << "Parsing file: " << filename << endl;
  ofstream myfile;
  myfile.open(filename, std::ios::out | std::ios::binary |std::ios::ate );
  if (myfile.is_open()){
    myfile.write((char*) &footprint,   sizeof(uint)*2);
    myfile.write((char*) &version,     sizeof(uint)*3);
    myfile.write((char*) &namelength,  sizeof(int));
    myfile.write(         filename,    sizeof(char) *namelength);
    for(int m=0;m<nmatrices;m++){
     // writeMatrixHeader(myfile,headers[m]);
      if(headers[m].type == INTMATRIX){ 
        writeintMatrixData(&myfile,headers[m],(*(intMatrix*)data[m]));
      }
      if(headers[m].type == DOUBLEMATRIX){
        writedoubleMatrixData(&myfile,headers[m],(*(doubleMatrix*)data[m]));
      }
      if(headers[m].type == FIXEDCHARMATRIX){
        writefixedCharMatrixData(&myfile,headers[m],(*(fixedcharMatrix*)data[m]));
      }
      if(headers[m].type == VARCHARMATRIX){
        writevarCharMatrixData(&myfile,headers[m],(*(varcharMatrix*)data[m]));
      }
      if(headers[m].type == MATRIXREFERENCE){
        writeReferenceMatrixData(&myfile,headers[m],data[m]);
      }
    }
    myfile.write((char*)&footprint, sizeof(uint)*2);
    myfile.close();
    return true;
  }else{
    Rprintf("File: %s could not be opened",filename);
    return false;
  }
  return false;
}
