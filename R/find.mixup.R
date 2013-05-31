#
# find.mixup.R
#
# copyright (c) 2012 Danny Arends
# last modified Jan, 2012
# first written Jan, 2012
# 

qtl.scan <- function(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), verbose = FALSE){
  if(missing(genotypes)) stop("genotypes are missing")
  if(missing(phenotypes)) stop("phenotypes are missing")
  results <- NULL
  rnames <- NULL
  cnt <- 1
  for(x in pheno.col){
    ss <- proc.time()
    results <- rbind(results,apply(genotypes,2, 
      function(geno){
        linmod <- lm(phenotypes[,x] ~ geno)
        -log10(anova(linmod)[[5]][1])
      }
    ))
    ee <- proc.time()
    if(verbose){
      cat("  - QTLscan of",colnames(phenotypes)[x],"took",as.numeric(ee[3]-ss[3]),"seconds\n")
    }
    rnames <- c(rnames,colnames(phenotypes)[x])
    cnt <- cnt +1
  }
  rownames(results) <- rnames
  class(results) <- c(class(results),"QTLscan")
  results
}

find.mixup <- function(phenotypes, genotypes, minLOD=15, nQTL=50, verbose = TRUE){
  if(missing(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)) stop("argument 'phenotypes' is missing, with no default")
  
  cat(" - Genotypes, ind=",nrow(genotypes),", markers=",ncol(genotypes),"\n",sep="")
  cat(" - Phenotypes, ind=",nrow(phenotypes),", phenotypes=",ncol(phenotypes),"\n",sep="")
  npheno <- ncol(phenotypes)
  nind   <- nrow(phenotypes)
  if(nQTL >= npheno) warning("too few traits in 'phenotypes'")
  scores <- rep(0,nind)
  it_cnt <- 1
  mQTL <- nQTL
  selected <- NULL
  used <- NULL
  while(nQTL > 0){
    toscan <- as.integer(runif(1)*npheno)
    while(toscan %in% used) toscan <- as.integer(runif(1)*npheno)
    qtlprofile <- qtl.scan(genotypes, phenotypes, toscan)
    if(max(abs(qtlprofile)) >= minLOD){
      marker <- which.max(abs(qtlprofile))
      nQTL <- (nQTL-1)
      cat(".")
      selected <- rbind(selected,c(toscan,marker))
    }else{
      #cat(toscan," ",max(abs(qtlprofile)),"\n")
    }
    used <- c(used,toscan)
    #cat("Itteration",it_cnt,", found",mQTL-nQTL,"QTLs",max(abs(qtlprofile)),"\n")
    it_cnt <- it_cnt+1
  }
  cat("\n")
  cat("Found enough QTL\n")
  mismatched <- NULL
  for(x in 1:nrow(selected)){
    g1 <- mean(phenotypes[genotypes[,selected[x,2]]==1,selected[x,1]],na.rm=T)
    g2 <- mean(phenotypes[genotypes[,selected[x,2]]==2,selected[x,1]],na.rm=T)
    cnt <- 1
    expected <- NULL
    for(ph in phenotypes[,selected[x,1]]){
      if(abs(ph-g1) < abs(ph-g2)){
        expected <- c(expected,1)
      }else{
        expected <- c(expected,2)
      }
    }
    mismatched <- rbind(mismatched,genotypes[,selected[x,2]]!=expected)
  }
  mismatched
}

toMixScore <- function(mismatched){
  mixscores <- apply(mismatched,2,sum,na.rm=T) / nrow(mismatched)
  mixscores
}

find.mixup.rqtl <- function(cross, minLOD=15, nQTL=50, verbose = TRUE){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  find.mixup(pull.pheno(cross),pull.geno(cross), minLOD, nQTL, verbose)
}

setwd("E:/GBIC/Konrad")
library(qtl)
load("DannyCross.rdata")
mu1 <- find.mixup.rqtl(cross,nQTL=100)
mu2 <- find.mixup.rqtl(cross,nQTL=100)
mu3 <- find.mixup.rqtl(cross,nQTL=100)
mu4 <- find.mixup.rqtl(cross,nQTL=100)
mu5 <- find.mixup.rqtl(cross,nQTL=100)
mu6 <- find.mixup.rqtl(cross,nQTL=100)
mu7 <- find.mixup.rqtl(cross,nQTL=100)
mu8 <- find.mixup.rqtl(cross,nQTL=100)
mu9 <- find.mixup.rqtl(cross,nQTL=100)
mu10 <- find.mixup.rqtl(cross,nQTL=100)

plot(toMixScore(rbind(mu1,mu2)),t='h')
