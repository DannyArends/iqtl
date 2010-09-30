#
# contrastqtlmapping.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: contrastqtlmapping, lodscorevectortoscanone
# Basic scripts for contrast QTL mapping using multiple regression
#

contrastqtlmapping <- function(cross,pheno.col=1,type=0,cofactors=NULL,verbose=TRUE){
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
  qtlprofile <- lodscorevectortoscanone(cross,lodscores)
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

lodscorevectortoscanone <- function(cross,lodscores){
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
  qtlprofile
}
