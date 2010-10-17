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
  result <- .C("inverseF_R",df1=as.integer(df1),
                            df2=as.integer(df2),
                            alpha=alpha,
                            out=0)
  result$out
}

backwardeliminate <- function(cross,pheno.col=1,cofactors,alpha=0.10){
  contrastlist  <- crosstocontrastlist(cross,0)
  contrastlist  <- lapply(FUN=scaledowncontrast,contrastlist)
  model         <- cofactors
  finished      <- FALSE
  pheno         <- cross$pheno[,pheno.col] 
  thismarkerclist <- contrastlisttodesignmatrix(contrastlist,model)
  logLfull      <- modellikelyhood(thismarkerclist,pheno)$likelihood
  dropneeded    <- 2*inverseF(2,nrow(thismarkerclist)-length(model),alpha);
  cat("Likelihood FULL model: ",logLfull,"\n");
  cat("Drop Needed: ",dropneeded,"\n");
  while(!finished && length(model) > 1){
    loglikelyhood <- NULL
    for(todrop in 1:length(model)){
      tempmodel <- model[-todrop]
      cat("Size temp model:",length(tempmodel),"/",length(model),"\n")
      thismarkerclist <- contrastlisttodesignmatrix(contrastlist,tempmodel)
      loglikelyhood = c(loglikelyhood,modellikelyhood(thismarkerclist,pheno)$likelihood)
    }
    leastinterestingmodel = which.max(loglikelyhood);
    likelihoodleastinterestingmodel = max(loglikelyhood);
    if(dropneeded > abs(logLfull - likelihoodleastinterestingmodel)){
      logLfull = likelihoodleastinterestingmodel;
      model <- model[-leastinterestingmodel]
      cat("Likelihood of the new full model: ",logLfull,"\n");
    }else{
      cat("\n\nWe have a model\n");
      for(x in 1:length(model)){
        cat("Variable ",model[x]," in Model\n");
      }
      finished=TRUE;
    }
  }
}