#
# backwardelimination.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: inverseF, backwardeliminate
#  Basic scripts for Causal inference
#

inverseF <- function(df1,df2,alpha=0.05){
  result <- .C("inverseF_R",df1=as.integer(df1),
                            df2=as.integer(df2),
                            alpha=alpha,
                            out=0)
  result$out
}

backwardeliminate <- function(cross,pheno.col=1,type=0,cofactors,alpha=0.01,verbose=FALSE){
  if(length(pheno.col) > 1){
    return(lapply(pheno.col,function(x){
      if(x <= nphe(cross) && x > 0){
        backwardeliminate(cross=cross,pheno.col=x,type=type,cofactors=cofactors,alpha=alpha,verbose=verbose)
      }
    }))
  }
  s <- proc.time()
  if(missing(cross)) stop("No cross object")
  if(missing(cofactors)) stop("No cofactors, please supply cofactors")
  contrastlist  <- crosstocontrastlist(cross,type)
  contrastlist  <- lapply(FUN=scaledowncontrast,contrastlist)
  out_of_range  <- which(cofactors > length(contrastlist) | cofactors < 1)
  if(length(out_of_range) > 0){
    warning("Some cofactors are out of range: ",cofactors[out_of_range])
    cofactors <- cofactors[-out_of_range];
  }
  model         <- cofactors
  finished      <- FALSE
  pheno         <- cross$pheno[,pheno.col]
  thismarkerclist <- contrastlisttodesignmatrix(contrastlist,model)
  logLfull      <- modellikelihood(thismarkerclist,pheno)$likelihood
  dropneeded    <- 2*inverseF(2,nrow(thismarkerclist)-length(model),alpha);
  e <- proc.time()
  startup <- as.numeric(e[3]-s[3])
  if(verbose){
    cat("Likelihood FULL model:",logLfull,"LR-cutoff:",dropneeded,"\n");
    cat("Full model in:",startup,"sec\n");
  }
  s <- proc.time()
  while(!finished && length(model) > 1){
    loglikelyhood <- NULL
    for(todrop in 1:length(model)){
      tempmodel <- model[-todrop]
      thismarkerclist <- contrastlisttodesignmatrix(contrastlist,tempmodel)
      loglikelyhood = c(loglikelyhood,modellikelihood(thismarkerclist,pheno)$likelihood)
    }
    leastinterestingmodel = which.max(loglikelyhood)
    likelihoodleastinterestingmodel = max(loglikelyhood)
    e <- proc.time()
    reduction <- as.numeric(e[3]-s[3])
    if(dropneeded > abs(logLfull - likelihoodleastinterestingmodel)){
      logLfull = likelihoodleastinterestingmodel;
      if(verbose) cat("Dropped variable",model[leastinterestingmodel],":",logLfull,"after",reduction,"sec\n");
      model <- model[-leastinterestingmodel]
    }else{
      finished=TRUE;
    }
  }
  if(!finished && dropneeded > abs(logLfull - min(loglikelyhood))){
    model <- model[-which.min(loglikelyhood)]
  }
  if(verbose){
    cat("Model: ");
    for(x in 1:length(model)){
      cat("Variable",model[x],"in Model\n");
    }
  }
  if(verbose){
    cat("Startup:",startup,"sec\n")
    cat("Modeling:",reduction,"sec\n")
  }
  model
}
