#
#
# read.xgap.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: read.xgap, XGAPheader, XGAPread
#
#


XGAPheader <- function(file){
	header <- .C("R_load_XGAPheader",
               filename=as.character(file),
               num_rows=as.integer(0),
               num_cols=as.integer(0),
               isNumeric=as.integer(0),
               pos=as.integer(0),
               PACKAGE="iqtl")
  if(header$num_rows < 0 || header$num_cols < 0) stop("No valid XGAP Header")
  result <- .C("R_load_XGAPnames",
               filename=as.character(file),
               num_rows=as.integer(header$num_rows),
               num_cols=as.integer(header$num_cols),
               pos=as.integer(header$pos),
               rownames=as.character(rep("",header$num_rows)),
               colnames=as.character(rep("",header$num_cols)),
               PACKAGE="iqtl")
  result$isNumeric = header$isNumeric
  result
}

XGAPread <- function(file, header){
  result <- NULL
  if(header$isNumeric){
    result <- .C("R_load_XGAPdouble",
                 filename=as.character(file),
                 num_rows=as.integer(header$num_rows),
                 num_cols=as.integer(header$num_cols),
                 pos=as.integer(header$pos),
                 data=as.numeric(matrix(0,header$num_rows,header$num_cols)),
                 PACKAGE="iqtl")
  }else{
    result <- .C("R_load_XGAPstring",
                 filename=as.character(file),
                 num_rows=as.integer(header$num_rows),
                 num_cols=as.integer(header$num_cols),
                 pos=as.integer(header$pos),
                 data=as.character(matrix(" ",header$num_rows,header$num_cols)),
                 PACKAGE="iqtl")
  }
  result$isNumeric=header$isNumeric
  result$rownames=header$rownames
  result$colnames=header$colnames
  result$data <- matrix(result$data,result$num_rows,result$num_cols,byrow=T)  
  result
}

read.xgap <- function(location=NULL){
  if(is.null(location)) stop("No such location")
  if(!(is.na(grep("http://",location)&&1))){
    content = getBinaryURL(location)
    tmp = tempfile()
    writeBin(content, con = tmp)
    location <- tmp
  }else{
    Tfile <- file(location, "r")
    close(Tfile)
  }
  h <- XGAPheader(location)
 
  d <- XGAPread(location,h)
  m <- d$data
  colnames(m) <- d$colnames
  rownames(m) <- d$rownames
  m
}
