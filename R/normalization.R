#
# normalization.R
#
# copyright (c) 2010, Danny Arends
# last modified Okt, 2010
# first written Okt, 2010
# 
# R functions: autonormalize
# Basic functions to do Normalization of the phenotypes in the cross object
# Depends on the library VGAM for probit
#

reciproce <- function(x){
  r <- 1/x
  r
}

donothing <- function(x){
  x
}

autonormalize <- function(cross,transformations=NULL,verbose=TRUE){
  originalfile <- cross
  newcross <- cross
  methods <- c(donothing,probit,log,sqrt,reciproce)
  chosen <- NULL
  if(verbose){
    cat("INFO: ------------------------------\n")
    cat("INFO: 1=No transformation\n")
    cat("INFO: 2=Probit transformation\n")
    cat("INFO: 3=Log transformation\n")
    cat("INFO: 4=sqrt\n")
    cat("INFO: 5=reciprocal\n")
    cat("INFO: ------------------------------\n")
  }
  for(x in 1:nphe(cross)){
    cnt <- 1
    suc <- FALSE
    if(!is.null(transformations) && length(transformations)==ncol(pull.pheno(cross)) && transformations[x]!=0){
      cross$pheno[[x]] <- methods[transformations[x]](originalfile$pheno[[x]])
      chosen <- c(chosen,transformations[x])
    }else{
      for(m in methods){
        if(!any(originalfile$pheno[[x]]<0,na.rm=T)){
          cross$pheno[[x]] <- m(originalfile$pheno[[x]])
        }else{
          #if(verbose) cat("[script] Method ",cnt," skipped for phenotype ",x,"\n")
        }
        if(any(is.na(cross$pheno[[x]])) || any(is.infinite(cross$pheno[[x]]))){
          #if(verbose) cat("[script] NA / Inf produced skipping ",cnt," for phenotype ",x,"\n")
        }else{
          if(var(cross$pheno[[x]])==0){
            #if(verbose) cat("[script] No variation, skipping ",cnt," for phenotype ",x,"\n")
          }else{
            if(!mqmtestnormal(cross, pheno.col = x)){
              #if(verbose) cat("[script] Phenotype ",x," NOT normal after ",cnt,"\n")
            }else{
              if(verbose) cat("Phenotype",x,"normal after",cnt,"\n")
              chosen <- c(chosen,cnt)
              newcross$pheno[[x]] <- cross$pheno[[x]]
              suc <- TRUE
              break;
            }
          }
        }
        cnt <- cnt+1
      }
      if(!suc){
        cat("Unable to normalize",x,"\n")
        chosen <- c(chosen,1)
        newcross$pheno[[x]] <- originalfile$pheno[[x]]
      }
    }
  }
  newcross$transformations <- chosen
  newcross
}
