
#testregression <-function(){
#  library(iqtl)
#  x <- matrix(runif(100),10,1)
#  w <- rep(1,10)
#  y <- x[,1] + runif(2)
#  o <- multipleregression(x,w,y)
#  o
#}

permarkercontrasts <- function(genotypes,m=1,verbose=FALSE){
  marker <- genotypes[,m]
  ncontrasts <- length(unique(marker))-1
  if(any(marker==666)){
    ncontrasts <- ncontrasts-1
  }
  values <- sort(unique(marker))
  contrastmatrix <- matrix(0,length(marker),ncontrasts)
  if(verbose) cat("Marker",m,"leads to: ",ncontrasts,"Contrasts\n")
  columnnames <- NULL
  for(x in 1:ncontrasts){
    base <- values[x]
    contrast <- values[x+1]
    for(m in 1:length(marker)){
      contrastmatrix[m,x] <- 0
      if(marker[m]==base){
        contrastmatrix[m,x] <- 1
      }
      if(marker[m]==contrast){
        contrastmatrix[m,x] <- -1
      }
    }
    columnnames <- c(columnnames,paste(base,":",contrast,sep=""))    
  }
  colnames(contrastmatrix) <- columnnames
  contrastmatrix
}

genomecontrasts <- function(genotypes,m=1,genomecontrasts,verbose=FALSE){
  marker <- genotypes[,m]
  ncontrasts <- length(unique(genomecontrasts))-1
  if(any(genomecontrasts==666)){
    ncontrasts <- ncontrasts-1
  }
  values <- sort(unique(genomecontrasts))
  contrastmatrix <- matrix(0,length(marker),ncontrasts)
  if(verbose) cat("Marker",m,"leads to: ",ncontrasts,"Contrasts\n")
  columnnames <- NULL
  for(x in 1:ncontrasts){
    base <- values[x]
    contrast <- values[x+1]
    for(m in 1:length(marker)){
      contrastmatrix[m,x] <- 0
      if(marker[m]==base){
        contrastmatrix[m,x] <- 1
      }
      if(marker[m]==contrast){
        contrastmatrix[m,x] <- -1
      }
    }
    columnnames <- c(columnnames,paste(base,":",contrast,sep=""))
  }
  colnames(contrastmatrix) <- columnnames
  contrastmatrix
}

crosstocontrastlist <- function(cross,type=0,verbose=FALSE){
  geno <- pull.geno(cross)
  if(any(is.na(geno))){
    if(verbose) cat("Missing values detected setting up fake values to capture it\n")
    geno[is.na(geno)] <- 666
  }
  res <- vector("list", ncol(geno))
  for(x in 1:ncol(geno)){
    res[[x]] <- 0
    if(type==0) res[[x]] <- as.matrix(genomecontrasts(geno,x,unique(sort(as.numeric(geno))),verbose))
    if(type==1) res[[x]] <- as.matrix(permarkercontrasts(geno,x,verbose))
  }
  res
}

scaledowncontrast <- function(contrastlistitem){
  res <- as.matrix(contrastlistitem[,apply(FUN=function(x) {length(unique(x))},contrastlistitem,MA=2)!=1])
  res
}

contrastlisttodesignmatrix <- function(contrastlist,m,cofactors = NULL,verbose=FALSE){
  designmatrix <- contrastlist[[m]]
  if(verbose) cat("Contrast matrix for marker",m,"\n")
  nullmatrixlayout <- rep(1,ncol(contrastlist[[m]]))
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


singlemarkercontrastmapping <- function(cross,pheno.col=1,type=0,cofactors=NULL,verbose=TRUE){
  s <- proc.time()
  contrastlist <- crosstocontrastlist(cross,type)
  contrastlist <- lapply(FUN=scaledowncontrast,contrastlist)
  pheno <- pull.pheno(cross)[,pheno.col]
  tokeep <- !is.na(pheno)
  pheno <- pheno[tokeep]
  lodscores <- NULL
  e <- proc.time()
  startup <- as.numeric(e[3]-s[3])
  prepare <- 0
  mapping <- 0
  for(x in 1:length(contrastlist)){
    s <- proc.time()
    thismarkerclist <- contrastlisttodesignmatrix(contrastlist,x,cofactors)
    contrastmarker <- as.matrix(thismarkerclist[[1]][tokeep,])
    nullmarkerlayout <- thismarkerclist[[2]]
    e <- proc.time()
    prepare <- prepare + as.numeric(e[3]-s[3])
    s <- proc.time()
    if(ncol(contrastmarker)>0){
      lodscores <- c(lodscores,multipleregression(contrastmarker,pheno,nullmodellayout=nullmarkerlayout)$likelihood)
    }else{
      lodscores <- c(lodscores,0)
    }
    e <- proc.time()
    mapping <- mapping + as.numeric(e[3]-s[3])
  }
  s <- proc.time()
  n <- unlist(lapply(FUN=colnames,pull.map(cross)))
  chr <- NULL
  if(!is.null(ncol(pull.map(cross)[[1]]))){
    d <- as.numeric(unlist(lapply(pull.map(cross),FUN=function(x) {x[1,]})))
    for(i in 1:nchr(cross)){chr <- c(chr,rep(names(cross$geno)[i], ncol(pull.map(cross)[[i]])))}
  }else{
    d <- as.numeric(unlist(pull.map(cross)))
    for(i in 1:nchr(cross)){chr <- c(chr,rep(names(cross$geno)[i], length(pull.map(cross)[[i]])))}
  }
  qtlprofile <- cbind(chr,d,lodscores)
  qtlprofile <- as.data.frame(qtlprofile)
  qtlprofile[,1] <- as.factor(chr)
  qtlprofile[,2] <- as.numeric(d)
  qtlprofile[,3] <- as.numeric(lodscores)
  rownames(qtlprofile) <- n
  class(qtlprofile) <- c("scanone", "data.frame")
  e <- proc.time()
  cleanup <- as.numeric(e[3]-s[3])
  if(verbose){
    cat("startup:",startup,"sec\n")
    cat("prepare:",prepare,"sec\n")
    cat("mapping:",mapping,"sec\n")
    cat("cleanup:",cleanup,"sec\n")
  }
  qtlprofile
}


multipleregression <- function(designmatrix,y,weight=rep(1,nrow(designmatrix)),nullmodellayout=rep(1,ncol(designmatrix)),verbose=FALSE){
  if(length(weight) != length(y)){
    stop("Not all samples have a weight, length(weight) != length(y) (Values:",length(weight)," != ",length(y),")")
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



