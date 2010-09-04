/**********************************************************************
 *
 * parser.h
 *
 * Copyright (c) 2010 Danny Arends
 *
 **********************************************************************/
#ifndef _PARSER_H
  #define _PARSER_H
  extern "C"
  {
  #include <R.h>
  #include "datastructures.h"
  
  typedef struct{
    uint            size;
    char*           memblock;
  }rawfilestruct;
  
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
