/**********************************************************************
 * 
 * xgapbin_r.cpp
 *
 * copyright (c) 2009 Danny Arends
 * last modified Sep, 2010
 * first written Apr, 2009
 *
 * C functions for the iqtl package
 * Contains: R_load_XGAPheader
 *           R_load_XGAPnames
 *           R_load_XGAPdouble
 *           R_load_XGAPstring
 *
 **********************************************************************/

extern "C"
{
#include "xgapbin_r.h"
#include "datastructures.h"
	
  void R_load_XGAPheader(char **filename,int* num_rows,int* num_cols,int* isNumeric,int* pos){
  //  Rprintf("File to load: %s\n",filename[0]);
    rawfilestruct binfile;
    if(parsefile(filename[0],&binfile)){
      char *matrixname,*investigationname,*colname,*rowname;
      uint newpos,nrows,ncols;
      uint deciortext;
      newpos = getname(binfile,&matrixname,1);
      newpos = getname(binfile,&investigationname,newpos);
      newpos = getname(binfile,&colname,newpos);
      newpos = getname(binfile,&rowname,newpos);
      deciortext = isbinary(binfile,newpos);
    //  Rprintf("Matrix= %s\n",matrixname);
    //  Rprintf("Numeric= %i\n",deciortext);
      newpos++;
      ncols = getncols(binfile,newpos);
    //  Rprintf("Cols= %i\n",ncols);
      newpos = newpos+4;
      nrows = getnrows(binfile,newpos);
    //  Rprintf("Rows= %i\n",nrows);
      *num_rows = nrows;
      *num_cols = ncols;
      *isNumeric = deciortext;
      newpos = newpos+4;

      (*pos) = newpos;      
      freememory(binfile);
    }
  }
  
void R_load_XGAPnames(char **filename,int* num_rows,int* num_cols,int* pos, char** rownames, char** colnames){
    rawfilestruct binfile;
    if(parsefile(filename[0],&binfile)){    
      uint newpos = (*pos);
      uint ncols = (*num_cols);
      uint nrows = (*num_rows);      
      uint *collengths,*rowlengths;
      collengths = getlengths(binfile,newpos,*num_cols);
      newpos = newpos+ (*num_cols);
      rowlengths = getlengths(binfile,newpos,*num_rows);
      newpos = newpos+ (*num_rows);
      char** cname;
      char** rname;
      newpos = newpos + getnames(&cname,binfile,newpos,*num_cols,collengths);
      newpos = newpos + getnames(&rname,binfile,newpos,*num_rows,rowlengths);
      //Rprintf("Parsed Flenames\n");
      for(uint r=0;r<nrows;r++){
        //Rprintf("%d->%s\n",r,rname[r]);
        rownames[r] = rname[r];
      }
      for(uint c=0;c<ncols;c++){
        colnames[c] = cname[c];
      }
      //Rprintf("In R mem\n");
      (*pos) = newpos;
      free(rname);
      free(cname);
      freememory(binfile);
    }
   
  }
  /*
      if(!deciortext){
        char*** data;      
        uint l = fixedstringlength(binfile,newpos);
        newpos++;
        if(l==0){
          //cout << "Variable lengths" << endl;
          uint *datalengths;
          datalengths = getlengths(binfile,newpos,nrows*ncols);
          getVariableStringdata(binfile,newpos,ncols,nrows,datalengths,data);
        }else{
         // cout << "Fixed length: " << l << endl;

          getFixedStringdata(binfile,newpos,ncols,nrows,l,data);
        }
      }else{
        double **data;
        getDecimaldata(binfile,newpos,ncols,nrows,data);
      }
      freememory(binfile);
    }
    */
    
  void R_load_XGAPdouble(char **filename,int* num_rows,int* num_cols,int* pos,double* data){
    rawfilestruct binfile;
    //Rprintf("Numeric File to load: %s\n",filename[0]);
    if(parsefile(filename[0],&binfile)){
      uint newpos = (*pos);
      uint ncols = (*num_cols);
      uint nrows = (*num_rows);
      double **rdata;
      getDecimaldata(binfile,newpos,ncols,nrows,&rdata);
      for(uint r=0;r<nrows;r++){
        for(uint c=0;c<ncols;c++){
          //Rprintf("%d.%d %d\n",r,c,(r*ncols)+c);
          data[(r*ncols)+c] = rdata[r][c];
          //Rprintf("Done: %d\n",(r*ncols)+c);
        }
      }
      for (uint i=0; i<nrows; i++) {
        free(rdata[i]);
      }
      free(rdata);
      freememory(binfile);      
      //Rprintf("Done ? \n");
    }
  }
  
void R_load_XGAPstring(char **filename,int* num_rows,int* num_cols,int* pos, char** data){
    rawfilestruct binfile;
  //  Rprintf("Text File to load: %s\n",filename[0]);
    if(parsefile(filename[0],&binfile)){
      uint newpos = (*pos);
     // Rprintf("pos:%d\n",newpos);    
      uint l = fixedstringlength(binfile,newpos);
      newpos++;
      uint ncols = (*num_cols);
     // Rprintf("col:%d\n",ncols);    
      uint nrows = (*num_rows);
      uint r=0;
      uint c=0;
     // Rprintf("row:%d\n",nrows);  

      char ***rdata;
      if(l > 0){
        //fixed string
        //Rprintf("Fixed\n");
        getFixedStringdata(binfile,newpos,ncols,nrows,l,&rdata);
      }else{
        //Variable String
        //Rprintf("Variable\n");
        uint* lengths;
        lengths = getlengths(binfile,newpos,ncols*nrows);
        newpos = newpos+(ncols*nrows);
        getVariableStringdata(binfile,newpos,ncols,nrows,lengths,&rdata);
      }
      for( r=0;r<nrows;r++){
        for( c=0;c<ncols;c++){
          //Rprintf("%d.%d %s\n",r,c,rdata[r][c]);
          data[(r*ncols)+c] = rdata[r][c];
        }
      }
      for( r=0;r<nrows;r++){
        free(rdata[r]);
      }
      Rprintf("Deallocated all rows \n");
      free(rdata);
      Rprintf("Deallocated matrix \n");
      freememory(binfile);      
      Rprintf("Free file struct \n");
      //Rprintf("Done ? \n");
    }
  }
}
