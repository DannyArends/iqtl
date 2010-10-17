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
  s <- proc.time()
  contrastlist  <- crosstocontrastlist(cross,type)
  contrastlist  <- lapply(FUN=scaledowncontrast,contrastlist)
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
    leastinterestingmodel = which.max(loglikelyhood);
    likelihoodleastinterestingmodel = max(loglikelyhood);
    e <- proc.time()
    reduction <- as.numeric(e[3]-s[3])
    if(dropneeded > abs(logLfull - likelihoodleastinterestingmodel)){
      logLfull = likelihoodleastinterestingmodel;
      if(verbose) cat("Dropped variable",model[leastinterestingmodel],":",logLfull,"after",reduction,"sec\n");
      model <- model[-leastinterestingmodel]
    }else{
      if(verbose){
        cat("We have a model\n");
        for(x in 1:length(model)){
          cat("Variable",model[x],"in Model\n");
        }
      }
      finished=TRUE;
    }
  }
  if(verbose){
    cat("startup:",startup,"sec\n")
    cat("modeling:",reduction,"sec\n")
  }
  model
}