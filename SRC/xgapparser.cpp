/**********************************************************************
 * 
 * xgapparser.cpp
 *
 * copyright (c) 2009 Danny Arends
 * last modified Sep, 2010
 * first written Apr, 2009
 *
 * C functions for the iqtl package
 * Contains:
 *
 **********************************************************************/
 
#include "xgapparser.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits> 
#include <cfloat>
#define LITTLE_ENDIAN 0
#define BIG_ENDIAN    1

int endian(void){
  //Checks endian setting of the machine
  int i = 1;
  char *p = (char *)&i;

  if (p[0] == 1){
    return LITTLE_ENDIAN;
  }else{
    return BIG_ENDIAN;
  }
}


char* integer2binary(uint a){
  //Used for debug prints integer as a binary
  char *str,*tmp;
  int cnt = 31;
  str = (char *) malloc(33); /*32 + 1 , becoz its a 32 bit bin number*/
  tmp = str;
  while ( cnt > -1 ){
    str[cnt]= '0';
    cnt --;
  }
  cnt = 31;
  while (a > 0){
    if (a%2==1){
      str[cnt] = '1';
    }
    cnt--;
    a = a/2 ;
  }
  return tmp;
}

void printcharasbinary(char ch){
  //Used for debug prints characters as a binary
  int i = CHAR_BIT;
  while (i > 0){
    -- i;
    //Rprintf("%d",(ch&(1 << i) ? '1' : '0'));
  }
}

bool parsefile(const char* filename, rawfilestruct* rawfile){
  //cout << "Parsing file: " << filename << endl;
  std::ifstream myfile;
  myfile.open(filename, std::ios::in | std::ios::binary |std::ios::ate );
  if (myfile.is_open()){
    rawfile->size = myfile.tellg();
    rawfile->memblock = (char*)mycalloc(rawfile->size,sizeof(char));
    myfile.seekg(0, std::ios::beg);
    myfile.read (rawfile->memblock, rawfile->size);
    myfile.close();
    //Rprintf("the complete file (%d Bytes) is in memory\n",rawfile->size);
    return true;
  }else{
    Rprintf("File: %s not found",filename);
    return false;
  }
  return false;
}

void rawwrite(rawfilestruct rawfile){
  for(uint x=0;x<=100;x++){
    printcharasbinary(rawfile.memblock[x]);
  }
}

char getnullchar(rawfilestruct rawfile){
  return rawfile.memblock[0];
}

bool isbinary(rawfilestruct rawfile, uint pos){
  return rawfile.memblock[pos];
}

uint getncols(rawfilestruct rawfile, uint pos){
  char* trans = (char*)mycalloc(3,sizeof(char));
  for(uint x=pos;x < pos+4;x++){
    trans[3-(x-pos)]= rawfile.memblock[x];
  }
  int t = *reinterpret_cast<int*>(trans);
  free(trans);
  return t;
}

uint getnrows(rawfilestruct rawfile, uint pos){
  char* trans = (char*)mycalloc(3,sizeof(char));
  for(uint x=pos;x < pos+4;x++){
    trans[3-(x-pos)]= rawfile.memblock[x];
  }
  int t = *reinterpret_cast<int*>(trans);
  free(trans);
  return t;
}

uint* getlengths(rawfilestruct rawfile, uint pos,uint n){
  uint* lengths = (uint*)mycalloc(n,sizeof(uint));
  for(uint x=0;x<n;x++){
    //int val = (int)rawfile.memblock[pos+x];
    lengths[x] = (uint)rawfile.memblock[pos+x];
    //No need we accepted 128 as maximum length
    //if(val < 0){
     // val = 128+ (128+val);
    //}
    //Rprintf("[%d] = %d -> %d\n",x,lengths[x],val);
  }
  //cout << endl;
  return lengths;
}

uint getnames(char*** names, rawfilestruct rawfile, uint pos,uint num, uint* sizes){
  (*names) = (char**)mycalloc((num),sizeof(char*));
  uint cnt = 0;
  for(uint n = 0; n < num;n++){
    char* tempname = (char*)mycalloc((sizes[n]+1),sizeof(char));
    for(uint t=0;t<sizes[n];t++){
        tempname[t] = rawfile.memblock[pos+cnt];
        cnt++;
    }
    tempname[sizes[n]] = '\0';
    (*names)[n] = tempname;
   // Rprintf("name %s\n",(*names)[n]);
    //cout << n << " name: " << names[n] << endl;
  }
  return cnt;
}

uint getname(rawfilestruct rawfile,char** name, uint pos){
  uint length = (uint)rawfile.memblock[pos];
  cout << "Name size: " << length << endl;
  pos++;
  *name = (char*)mycalloc(length+1,sizeof(char));
  for(uint i=0; i<length; i++){
    (*name)[i] = rawfile.memblock[pos+i];
  }
  (*name)[length] = '\0';
  //cout << "Found: " << length+pos << " " << name << endl;
  return (pos+length);
}

uint fixedstringlength(rawfilestruct rawfile,uint pos){
  uint length;
  length = rawfile.memblock[pos];
  return length;
}

void getDecimaldata(rawfilestruct rawfile,uint pos,uint cols,uint rows, double ***data){

  (*data) = (double**)mycalloc(rows,sizeof(double*));
  for(uint i=0;i<rows;i++){
    (*data)[i] = (double*)mycalloc(cols,sizeof(double));
  }
  for(uint r=0;r<rows;r++){    
    for(uint c=0;c<cols;c++){
      char* trans = (char*)mycalloc(8,sizeof(char));      
      for(uint x=pos;x < pos+8;x++){
        trans[7-(x-pos)]= rawfile.memblock[x];
      }
      (*data)[r][c] = *reinterpret_cast<double*>(trans);
      if((*data)[r][c] >= std::numeric_limits<double>::max() ){
        (*data)[r][c] = NAN;
     //   Rprintf("NaN\t");
      }else{
    //    Rprintf("%.3f\t",(*data)[r][c]);
      }
      free(trans);
      pos = pos+8;
    }
  //  Rprintf("\n");
  }
}


void getFixedStringdata(rawfilestruct rawfile,uint pos,uint cols,uint rows,uint l,char ****data){
//  Rprintf("FixStringData\n");
  (*data) = (char***)mycalloc(rows,sizeof(char**));
  for(uint i=0;i<rows;i++){
    (*data)[i] = (char**)mycalloc(cols,sizeof(char*));
  }
  for(uint r=0;r<rows;r++){    
    for(uint c=0;c<cols;c++){
     // Rprintf("%d %d\n",r,c);    
      char *n = (char*)mycalloc(l+1,sizeof(char));
      for(uint i=0; i<l; i++){
        n[i] = rawfile.memblock[pos+i];
      }
      n[l] = '\0';
    //  Rprintf("%d %d %s\n",r,c,n);
      (*data)[r][c] = n;
      pos=pos+l;
      free(n);
    }
  }
  //Rprintf("Done\n");
}

void getVariableStringdata(rawfilestruct rawfile,uint pos,uint cols,uint rows,uint* l, char ****data){
  //Rprintf("VarStringData\n");
  (*data) = (char***)mycalloc(rows,sizeof(char**));
  for(uint i=0;i<rows;i++){
    (*data)[i] = (char**)mycalloc(cols,sizeof(char*));
  }
  uint cnt=0;
  for(uint r=0;r<rows;r++){    
    for(uint c=0;c<cols;c++){
      //Rprintf("%d %d-> %d\n",r,c,l[cnt]+1);
      char *n = (char*)mycalloc(l[cnt]+1,sizeof(char));
      for(uint i=0; i<l[cnt]; i++){
        //Rprintf("%d\n",i); 
        n[i] = rawfile.memblock[pos+i];
      }
      n[l[cnt]] = '\0';
      //Rprintf("%d %d %s\n",r,c,n);
      (*data)[r][c] = n;
      cnt++;
      free(n);
      pos=pos+l[cnt];
    }
  }
 // Rprintf("Done\n");
}


void freememory(rawfilestruct rawfile){
  rawfile.size= 0;
  free(rawfile.memblock);
}

