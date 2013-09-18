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
  no_y <- which(is.na(y))
  if(length(no_y)>0){
    if(verbose) cat("Dropping individuals (",paste(which(is.na(y))),")with missing response values\n")
    y <- y[-no_y]
    weight <- weight[-no_y]
    designmatrix <- as.matrix(designmatrix[-no_y,],dim(designmatrix)[0],dim(designmatrix)[1])
  }
  if(nrow(designmatrix) != length(weight)){
    stop("Not all samples have a weight, nrow(designmatrix) != length(weight) (Values:",length(weight)," != ",length(y),")")
  }
  if(nrow(designmatrix) != length(y)){
    stop("Not all samples have a output, nrow(designmatrix) != length(y) (Values:",nrow(designmatrix)," != ",length(y),")")
  }
  if(!sum(designmatrix[,1])==length(y)){
    if(verbose) warning("Adding estimate of constant in model")
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

modellikelihood <- function(designmatrix,y,weight=rep(1,nrow(designmatrix)),verbose=FALSE){
  no_y <- which(is.na(y))
  if(length(no_y)>0){
    if(verbose) cat("Dropping individuals (",paste(which(is.na(y))),")with missing response values\n")
    y <- y[-no_y]
    weight <- weight[-no_y]
    designmatrix <- as.matrix(designmatrix[-no_y,],dim(designmatrix)[0],dim(designmatrix)[1])
  }
  if(nrow(designmatrix) != length(weight)){
    stop("Not all samples have a weight, nrow(designmatrix) != length(weight) (Values:",nrow(designmatrix)," != ",length(weight),")")
  }
  if(nrow(designmatrix) != length(y)){
    stop("Not all samples have a output, nrow(designmatrix) != length(y) (Values:",nrow(designmatrix)," != ",length(y),")")
  }
  if(!sum(designmatrix[,1])==length(y)){
    if(verbose) warning("Adding estimate of constant in model")
    designmatrix <- cbind(rep(1,length(y)),designmatrix)
  }
  result <- .C("modellikelihoodbyem_R",nvariables=as.integer(ncol(as.matrix(designmatrix))),
                                nsamples=as.integer(nrow(as.matrix(designmatrix))),
                                x=as.matrix(designmatrix),
                                w=weight,
                                y=y,
                                verbose=as.integer(verbose),
                                likelihood=0)
  result
}
