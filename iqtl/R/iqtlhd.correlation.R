#
# iqtlhd.correlation.R
#
# copyright (c) 2012 Danny Arends
# last modified Jan, 2012
# first written Oct, 2011
# 

correlation <- function(x, y, verbose = FALSE){
  if(is.matrix(x) && dim(x)[2] > 1){
    return(matrix_correlation(x,verbose=verbose))
  }else{
    if(!has_qtlHD()){
      if(verbose) cat(.has_qtlHD_warnmsg)
      return(cor(x,y))
    }else{
      if(missing(y))stop("'y' is missing");
      if(length(x) == length(y)){
        return(.C("correlation_v",
          x=as.numeric(x),
          y=as.numeric(y),
          as.integer(length(x)),
          cor=as.numeric(0),
          PACKAGE="qtlHD")$cor
        )
      }else{
        stop("The length of 'x' is not equal to the length of 'y'")
      }
    }
  }
}

matrix_correlation <- function(x, verbose = FALSE){
  if(!has_qtlHD()){
    if(verbose) cat(.has_qtlHD_warnmsg)
    return(cor(x))
  }else{
      res <- .C("correlation_m",
      nrx=as.integer(nrow(x)),
      ncx=as.integer(ncol(x)),
      x=as.numeric(x),
      cor=as.numeric(rep(1,ncol(x)*ncol(x))),
      PACKAGE="qtlHD"
    )
    mcor = matrix(res$cor,ncol(x),ncol(x))
    return(mcor)
  }
}
