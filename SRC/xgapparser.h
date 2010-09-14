/**
 * \file xgapparser.h
 * \brief Header file, for xgapparser.cpp
 *
 * last modified Sep, 2010
 * first written Apr, 2009
 * copyright (c) 2009-2010 Danny Arends
 *
 * C functions for the iqtl package
 * Contains:  endian
 *            integer2binary
 *            printcharasbinary
 *            parsefile
 *            getnullchar
 *            isbinary
 *            getncols
 *            getnrows
 *
 **/

#ifndef _XGAPPARSER_H
  #define _XGAPPARSER_H
  extern "C"
  {
  #include <R.h>
  #include "datastructures.h"
  
  typedef struct{
    uint            size;
    char*           memblock;
  }rawfilestruct;
  
  int endian(void);
  char* integer2binary(uint a);
  void printcharasbinary(char ch);
  bool parsefile(const char* filename, rawfilestruct* rawfile);
  char getnullchar(rawfilestruct rawfile);
  bool isbinary(rawfilestruct rawfile, uint pos);
  uint getncols(rawfilestruct rawfile, uint pos);
  uint getnrows(rawfilestruct rawfile, uint pos);
  uint* getlengths(rawfilestruct rawfile, uint pos,uint n);
  uint getnames(char*** names, rawfilestruct rawfile, uint pos,uint num, uint* sizes);
  uint getname(rawfilestruct rawfile,char** name, uint pos);
  uint fixedstringlength(rawfilestruct rawfile,uint pos);
  void getDecimaldata(rawfilestruct rawfile,uint pos,uint cols,uint rows, double ***data);
  void getFixedStringdata(rawfilestruct rawfile,uint pos,uint cols,uint rows,uint l,char ****data);
  void getVariableStringdata(rawfilestruct rawfile,uint pos,uint cols,uint rows,uint* l, char ****data);
  void rawwrite(rawfilestruct rawfile);
  void freememory( rawfilestruct rawfile);
  
  }
#endif
