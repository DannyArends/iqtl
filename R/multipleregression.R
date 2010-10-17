#
# multipleregression.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: multipleregression
# Basic scripts for Multiple Regression 
#
#testregression <-function(){
#  library(iqtl)
#  x <- matrix(runif(100),10,1)
#  w <- rep(1,10)
#  y <- x[,1] + runif(2)
#  o <- multipleregression(x,y,w)
#  o
#}

multipleregression <- function(designmatrix,y,weight=rep(1,nrow(designmatrix)),nullmodellayout=rep(1,ncol(designmatrix)),verbose=FALSE){
  if(nrow(designmatrix) != length(weight)){
    stop("Not all samples have a weight, nrow(designmatrix) != length(weight) (Values:",length(weight)," != ",length(y),")")
  }
  if(nrow(designmatrix) != length(y)){
    stop("Not all samples have a output, nrow(designmatrix) != length(y) (Values:",nrow(designmatrix)," != ",length(y),")")
  }
  if(!sum(designmatrix[,1])==length(y)){
    warning("Adding estimate of constant in model")
    designmatrix <- cbind(rep(1,length(y)),designmatrix)
  }
  result <- .C("lodscorebyem_R",nvariables=as.integer(ncol(as.matrix(designmatrix))),
                                nsamples=as.integer(nrow(as.matrix(designmatrix))),
                                x=as.matrix(designmatrix),
                                w=weight,
                                y=y,
                                estparams=rep(1.0,as.integer(ncol(as.matrix(designmatrix)))),
                                nullmodellayout=as.integer(nullmodellayout),
                                verbose=as.integer(verbose),
                                likelihood=0)
  result
}

backwardelimination <- function(designmatrix,weight,y,verbose=FALSE){
  if(length(weight) != length(y)){
    stop("Not all samples have a weight, length(weight) != length(y)")
  }
  if(nrow(designmatrix) != length(y)){
    stop("Not all samples have a output, nrow(designmatrix) != length(y)")
  }
  if(!sum(designmatrix[,1])==length(y)){
    warning("Adding estimate of constant in model")
    designmatrix <- cbind(rep(1,length(y)),designmatrix)
  }
  result <- .C("backwardelimination_R",nvariables=as.integer(ncol(as.matrix(designmatrix))),
                                  nsamples=as.integer(nrow(as.matrix(designmatrix))),
                                  x=as.matrix(designmatrix),
                                  w=weight,
                                  y=y,
                                  verbose=as.integer(verbose),
                                  likelihood=0)
  result
}

#for(x in 1:2){testregression()}



