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
 
#include "xgapparser.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits> 
#include <cfloat>
#define LITTLE_ENDIAN 0
#define BIG_ENDIAN    1


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

void writeMatrixHeader(ofstream myfile, matrixHeader h){
  myfile.write((char*) &h.namelength,  sizeof(int));
  myfile.write(         h.name,        sizeof(char) * h.namelength);
  myfile.write((char*) &h.type,        sizeof(DataType));
  myfile.write((char*) &h.columns,     sizeof(int));
  myfile.write((char*) &h.rows,        sizeof(int));
}

void writeNumericalMatrixData(ofstream myfile, matrixHeader h, void* data){
  if(h.type == INTMATRIX){
    intMatrix* matrix = (intMatrix*)data;
  }
  if(h.type == DOUBLEMATRIX){
    doubleMatrix* matrix = ((doubleMatrix*)data);
  }
  for(int r;r<h.rows;r++){
    for(int c;c<h.columns;c++){  
        myfile.write((char*) &(*matrix).data[r][c], sizeof(matrix[r][c]));
      }
    }
  }
}

void writeCharMatrixData(ofstream myfile, matrixHeader h, void* data){
  if(h.type == FIXEDCHARMATRIX){
    fixedcharMatrix* matrix = (fixedcharMatrix*)data;
    myfile.write((char*) &(*matrix).length, sizeof(int));
  }
  if(h.type == VARCHARMATRIX){
    varcharMatrix* matrix = ((varcharMatrix*)data;
    for(int indx=0;indx < ;indx++){
      myfile.write((char*) &(*matrix).lengths, h.rows*h.cols* sizeof(int));
  }
  for(int r;r<h.rows;r++){
    for(int c;c<h.columns;c++){  
        myfile.write((char*) &(*matrix).data[r][c], sizeof(matrix[r][c]));
      }
    }
  }
}

void writeReferenceMatrixData(ofstream myfile, matrixHeader h, void* data){

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
      writeMatrixHeader(myfile,headers[m]);
      if(headers[m].type == INTMATRIX || headers[m].type == DOUBLEMATRIX){
        writeNumericalMatrixData(myfile,headers[m],data[m]);
      }
      if(headers[m].type == FIXEDCHARMATRIX|| headers[m].type == VARCHARMATRIX){
        writeCharMatrixData(myfile,headers[m],data[m]);
      }
      if(headers[m].type == MATRIXREFERENCE){
        writeReferenceMatrixData(myfile,headers[m],data[m]);
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