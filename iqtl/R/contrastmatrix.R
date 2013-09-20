#
# contrastmatrix.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: crosstocontrastlist,scaledowncontrast, contrastlisttodesignmatrix, permarkercontrasts, genomecontrasts
# Basic scripts for supporting contrast QTL mapping using multiple regression
#

crosstocontrastlist <- function(cross, type = 0, verbose = FALSE){
  geno <- pull.geno(cross)
  if(any(is.na(geno))){
    if(verbose) cat("Missing values detected setting up fake values to capture it\n")
    geno[is.na(geno)] <- 666
  }
  res <- vector("list", ncol(geno))
  for(x in 1:ncol(geno)){
    res[[x]] <- 0
    if(type==0) res[[x]] <- as.matrix(genomecontrasts(geno,x,verbose))
    if(type==1) res[[x]] <- as.matrix(markercontrasts(geno,x,verbose))
  }
  res
}

getThird <- function(x){ return(x[,3]) }

scaledowncontrast <- function(contrastlistitem){
  return( as.matrix(contrastlistitem[,apply(contrastlistitem, MARGIN = 2, FUN=function(x){ 
  	length(unique(x)) }) != 1] ) 
  )
}

contrastlisttomatrices <- function(contrastlist,m,cofactors = NULL,verbose=FALSE){
  designmatrix <- contrastlist[[m]]
  if(verbose) cat("Contrast matrix for marker",m,"\n")
  nullmatrixlayout <- rep(1, ncol(contrastlist[[m]]))
  if(!is.null(cofactors)){
    for(x in cofactors){
      #COM: Should be inrange function (or correlated)
      if(abs(cor(contrastlist[[x]],contrastlist[[m]])) < 0.6 ){
        designmatrix <- cbind(contrastlist[[x]],designmatrix)
        nullmatrixlayout <- c(rep(0,ncol(contrastlist[[x]])),nullmatrixlayout)
      }
    }
  }
  list(designmatrix,nullmatrixlayout)
}

contrastlisttodesignmatrix <- function(contrastlist, cofactors,verbose=FALSE){
  if(missing(contrastlist)) stop("No contrastlist, please create one using crosstocontrastlist")
  if(missing(cofactors)) stop("No cofactors, please supply cofactors")
  designmatrix <- NULL
  if(!is.null(cofactors)){
    for(x in cofactors){
      designmatrix <- cbind(contrastlist[[x]],designmatrix)
    }
  }
  designmatrix
}

markercontrasts <- function(genotypes,m=1,verbose=FALSE){
  marker <- genotypes[,m]
  ncontrasts <- length(unique(marker))-1
  if(any(marker==666)){
    ncontrasts <- ncontrasts-1
  }
  values <- sort(unique(marker))
  contrastmatrix <- matrix(0,length(marker),ncontrasts)
  if(verbose) cat("Marker", m, ":", ncontrasts, "contrasts\n")
  columnnames <- NULL
  for(x in 1:ncontrasts){
    base <- values[x]
    contrast <- values[x+1]
    for(m in 1:length(marker)){
      contrastmatrix[m,x] <- 0
      if(marker[m]==base){ contrastmatrix[m,x] <- 1 }
      if(marker[m]==contrast){ contrastmatrix[m,x] <- -1 }
    }
    columnnames <- c(columnnames,paste(base, ":",contrast,sep=""))    
  }
  colnames(contrastmatrix) <- columnnames
  contrastmatrix
}

genomecontrasts <- function(genotypes,m=1,verbose=FALSE){
  genomecontrasts <- unique(sort(as.numeric(genotypes)))
  ncontrasts <- length(unique(genomecontrasts))-1
  marker <- genotypes[,m]
  if(any(genomecontrasts==666)){ ncontrasts <- ncontrasts-1 }
  values <- sort(unique(genomecontrasts))
  contrastmatrix <- matrix(0,length(marker),ncontrasts)
  if(verbose) cat("Marker", m, ":", ncontrasts, "contrasts\n")
  columnnames <- NULL
  for(x in 1:ncontrasts){
    base <- values[x]
    contrast <- values[x+1]
    for(m in 1:length(marker)){
      contrastmatrix[m,x] <- 0
      if(marker[m]==base){ contrastmatrix[m,x] <- 1 }
      if(marker[m]==contrast){ contrastmatrix[m,x] <- -1 }
    }
    columnnames <- c(columnnames,paste(base,":",contrast,sep=""))
  }
  colnames(contrastmatrix) <- columnnames
  contrastmatrix
}
