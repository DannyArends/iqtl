#
# contrastqtlmapping.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: contrastqtlmapping, contrastqtlsignificance, getlodthreshold, lodscorestoscanone,
# Basic scripts for contrast QTL mapping using multiple regression
#

#Function used for estimating significance
contrastqtlsignificance <- function(cross, pheno.col = 1, type = 0, cofactors=NULL, cycles = 1000, verbose=TRUE){
  s <- proc.time()
  samples <- nrow(cross$pheno)
  originalphenotype <- cross$pheno[,pheno.col]
  lodscorematrix <- NULL
  #Pre-compute the contrast list to use the contrastqtlmapping.internal
  contrastlist <- crosstocontrastlist(cross,type)
  contrastlist <- lapply(FUN=scaledowncontrast,contrastlist)
  for(x in 1:cycles){
    cross$pheno[,pheno.col] <- cross$pheno[sample(samples),pheno.col]
    qtlresults <- contrastqtlmapping.internal(cross,contrastlist,pheno.col=pheno.col,type=type,cofactors=cofactors,verbose=FALSE)
    lodscorematrix <- rbind(lodscorematrix, qtlresults$lod)
  }
  cross$pheno[,pheno.col] <- originalphenotype
  e <- proc.time()
  if(verbose){
    cat("Total time:", as.numeric(e[3]-s[3]),"sec\n")
  }
  lodscorematrix
}

#Main function used for testing the new multiple regression interface to R/qtl
contrastqtlmapping <- function(cross,pheno.col=1,type=0,cofactors=NULL,verbose=TRUE){
  contrastlist <- crosstocontrastlist(cross,type)
  contrastlist <- lapply(FUN=scaledowncontrast,contrastlist)
  contrastqtlmapping.internal(cross,contrastlist,pheno.col,type,cofactors,verbose)
}

#Pre-compute the contrast list to use the contrastqtlmapping.internal, 
#Speeds up computation when doing repeated mappings on the same cross object
#
#Example:
# library(iqtl)
# data(hyper)
# contrastlist <- crosstocontrastlist(hyper,0)
# contrastlist <- lapply(FUN=scaledowncontrast,contrastlist)
# #Repeated internal call
# contrastqtlmapping.internal(hyper,contrastlist,pheno.col=1:1000)
  
contrastqtlmapping.internal <- function(cross,contrastlist=crosstocontrastlist(cross,type),pheno.col=1,type=0,cofactors=NULL,verbose=TRUE){
  s <- proc.time()
  phenotypes <- as.matrix(pull.pheno(cross)[,pheno.col],nind(cross),length(pheno.col))
  lodmatrix <- NULL
  e <- proc.time()
  startup <- as.numeric(e[3]-s[3])
  prepare <- 0
  mapping <- 0
  for(x in 1:length(contrastlist)){
    s <- proc.time()
    thismarkerclist <- contrastlisttomatrices(contrastlist,x,cofactors)
    nullmarkerlayout <- thismarkerclist[[2]]
    lodscores <- NULL
    e <- proc.time()
    prepare <- prepare + as.numeric(e[3]-s[3])
    
    for(p in 1:ncol(phenotypes)){
      s <- proc.time()
      tokeep <- !is.na(phenotypes[,p])
      
      pheno <- phenotypes[tokeep,p]
      contrastmarker <- as.matrix(thismarkerclist[[1]][tokeep,])
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
    lodmatrix <- rbind(lodmatrix,lodscores)
  }
  s <- proc.time()
  qtlprofile <- lodscorestoscanone(cross,lodmatrix)
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

#Gets the lod threshold from a matrix where each row is a permutation scan
plotlodscorematrix <- function(cross, result, lodscorematrix){
  plot(result,lodscorestoscanone(cross,rep(getlodthreshold(lodscorematrix,10),sum(nmar(cross)))),lodscorestoscanone(cross,rep(getlodthreshold(lodscorematrix),sum(nmar(cross)))),col=c("black","yellow","green"))
  legend("topright",c("Trait","5%","10%"),col=c("black","green","yellow"),lwd=c(1,1,1))
}

#Gets the lod threshold from a matrix where each row is a permutation scan
#unprecise when there are few lodscores ( < 500 )
getlodthreshold <- function(lodscorematrix,percentage = 5){
  if(!percentage > 0 && !percentage < 100){
    stop("Unknown percentage, must be 0 < percentage < 100")
  }
  sort(apply(lodscorematrix,1,max))[floor(nrow(lodscorematrix)*((100-percentage)/100))]
}

#To use karl's pretty histogram plot function
lodscorematrixtoscanoneperm <- function(lodscorematrix){
  lodscorematrix <- as.matrix(lodscorematrix)
  colnames(lodscorematrix) <- "lod"
  rownames(lodscorematrix) <- 1:ncol(lodscorematrix)
  class(lodscorematrix) <- c("scanoneperm",class(lodscorematrix))
  lodscorematrix
}

#Change any list of lodscores into a scanone object (only pre-req: length(lodscores)==sum(nmar(cross))
lodscorestoscanone <- function(cross,lodscores,traitnames = NULL){
  n <- unlist(lapply(FUN=names,pull.map(cross)))
  chr <- NULL
  if(!is.null(ncol(pull.map(cross)[[1]]))){
    d <- as.numeric(unlist(lapply(pull.map(cross),FUN=function(x) {x[1,]})))
    for(i in 1:nchr(cross)){
      chr <- c(chr,rep(names(cross$geno)[i], ncol(pull.map(cross)[[i]])))
    }
  }else{
    d <- as.numeric(unlist(pull.map(cross)))
    for(i in 1:nchr(cross)){
      chr <- c(chr,rep(names(cross$geno)[i], length(pull.map(cross)[[i]])))
    }
  }
  qtlprofile <- cbind(chr,d,lodscores)
  qtlprofile <- as.data.frame(qtlprofile)
  qtlprofile[,1] <- chr
  qtlprofile[,2] <- as.numeric(d)
  if(!is.null(ncol(lodscores))){
    for(x in 1:ncol(lodscores)){
      qtlprofile[,2+x] <- as.numeric(lodscores[,x])
    }
    traitnames = Names("lod",ncol(lodscores))
  }else{
     qtlprofile[,3] <- as.numeric(lodscores)
     traitnames = "lod"
  }
  rownames(qtlprofile) <- n
  colnames(qtlprofile) <- c("chr","cM",traitnames)
  class(qtlprofile) <- c("scanone", "data.frame")
  qtlprofile
}
