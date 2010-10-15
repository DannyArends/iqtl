#
# backwardelimination.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: mqmmultitoscanone, mqmmodelsasnetwork
#  Basic scripts for Causal inference
#

inverseF <- function(df1,df2,alpha=0.05){
  m <- NULL
  for(x in 1:length(df1)){
    r <- NULL
    for(y in 1:length(df2)){
    result <- .C("inverseF_R",df1=as.integer(df1[x]),
                              df2=as.integer(df2[y]),
                              alpha=alpha,
                              out=0)
    r <- c(r,result$out)
    }
    m <- rbind(m,r)
  }
  rownames(m) <- df1
  colnames(m) <- df2
  m
}

backwardeliminate <- function(cross,cofactors,alpha=0.05){
  contrastlist  <- crosstocontrastlist(cross,0)
  cat("1\n")
  contrastlist  <- lapply(FUN=scaledowncontrast,contrastlist)
  cat("2\n")
  model         <- cofactors
  cat("3\n")
  thismarkerclist <- contrastlisttodesignmatrix(contrastlist,cofactors)
  cat("4\n")
  logLfull      <- modellikelyhood(thismarkerclist, model);
  cat("5\n")
  dropneeded    <- 2*inverseF(2,nsamples-nvariables,alpha);
  while(!finished && sum(model) > 1){
    loglikelyhood <- vector("list", sum(model))
    for(todrop in 1:sum(model)){
      tempmodel <- model[-todrop]
      thismarkerclist <- contrastlisttodesignmatrix(contrastlist,cofactors)
      loglikelyhood[todrop] = modellikelyhood(thismarkerclist, model);
    }
    leastinterestingmodel = which.max(loglikelyhood);
    likelihoodleastinterestingmodel = max(loglikelyhood);
    if(dropneeded > abs(logLfull - likelihoodleastinterestingmodel)){
      logLfull = likelihoodleastinterestingmodel;
      if(verbose) cat("Likelihood of the new full model: ",logLfull,"\n");
    }else{
      cat("\n\nWe have a model\n");
      for(x in 1:sum(cofactors)){
        if(model[x]) cat("Variable ",x," in Model\n");
      }
      finished=true;
    }
  }
}