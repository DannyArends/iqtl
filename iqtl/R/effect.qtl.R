# QTLeffectsCorrection
# Removal of upto 5 main QTL effects in a RIL crossobject
# (c) 2010 Danny Arends

markereffect <- function(markergenotype,trait){
  class1 <- mean(trait[which(markergenotype==1)], na.rm=T)
  class2 <- mean(trait[which(markergenotype==2)], na.rm=T)
  oamean <- mean(trait, na.rm=T)
  return(c(class1-oamean, class2-oamean))
}

markereffects <- function(genotypes,trait){ apply(genotypes, 2, markereffect, trait=trait) }

individualcorrection <- function(trait, genotypes, verbose=FALSE){
  otrait <- trait
  eff <- markereffects(genotypes,trait)
  pprofile <- apply(genotypes,2,function(x){as.numeric(anova(lm(trait ~ x))[[5]][1])})
  indeffects <- rep(0,nrow(genotypes))
  prevmin <- -5
  maxqtl <- 5
  n <- -1
  while(min(pprofile,na.rm=T) < 0.001 && n != which.min(pprofile) && maxqtl > 0){
    n <- which.min(pprofile)
    indeffects <- indeffects + eff[genotypes[,n],n]
    trait <- otrait-indeffects
    tryCatch(
      pprofile <- apply(genotypes,2,function(x){
      	as.numeric(anova(lm(trait ~ x,na.action=na.exclude))[[5]][1])
      })
      ,error = function(x){cat(x[[1]],"\n");return(indeffects)}
    )
    maxqtl<- maxqtl - 1
  }
  if(verbose){
    if(maxqtl <= 0){
      cat("Removed 5 QTLs")
    }else{
      cat("Removed ", 5-maxqtl, "QTLs")
    }
  }
  indeffects
}

correcttrait <- function(cross, pheno.col=1, plotprofiles=TRUE, verbose=FALSE){
  if(plotprofiles) old <- scanone(cross,pheno.col=pheno.col)
  correction <- individualcorrection(trait=as.numeric(as.character(pull.pheno(cross)[,pheno.col])),genotypes=pull.geno(cross),verbose=verbose)
  if(plotprofiles){
    cross$pheno[,pheno.col] <- cross$pheno[,pheno.col]-correction
    new <- scanone(cross,pheno.col=pheno.col)
    plot(old,new)
  }
  correction
}

#SNOW enabled function to reduce computational burden
correctQTLeffect <- function(cross, pheno.col=NULL, n.clusters=3, batchsize=30, plotprofiles=FALSE, verbose=TRUE){
  totpheno <- pheno.col
  if(is.null(pheno.col)){ totpheno <- 1:ncol(cross$pheno) }
  cross <- calc.genoprob(fill.geno(cross))
  l <- as.list(totpheno)
  res <- NULL
  if(n.clusters > 1 && verbose) cat("-- Using SNOW package for computation\n")
  for(x in seq(1,length(totpheno),batchsize)){
    s <- proc.time()
    cat("Doing:",x,"..",min(x+(batchsize-1),length(totpheno)),"\n")
    if(n.clusters > 1){
      cl <- makeCluster(rep("localhost",n.clusters))
      r <- parLapply(cl,l[x:min(x+(batchsize-1),length(totpheno))],function(x){library(iqtl);correcttrait(cross,x,plotprofiles)})
      stopCluster(cl)
    }else{
      r <- lapply(l[x:min(x+(batchsize-1),length(totpheno))],function(x){library(iqtl);correcttrait(cross,x,plotprofiles)})
    }
    res <- c(res,r)
    e <- proc.time()
    if(batchsize <= length(pheno.col)){
      cat(batchsize, "phenotypes took:", as.numeric(e[3]-s[3]), "secs\n")
    }else{
      cat(length(totpheno), "phenotypes took:",as.numeric(e[3]-s[3]), "secs\n")
    }
  }
  invisible(matrix(unlist(res),length(res[[1]]),length(res),byrow=F))
}
