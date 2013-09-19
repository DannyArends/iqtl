#
#
# traitCorrelation.R
#
# copyright (c) 2011 Danny Arends
# last modified apr, 2011
# first written apr, 2011
# 
# R functions to calculate correlations between 2 cross objects (Assumption: Holding the same individuals)
#
# Contains: traitCorrelation
#

compareTraits <- function(cross1, cross2, FUN = cor, verbose=TRUE){
  if(nind(cross1) != nind(cross2)) stop("No match between number of individuals")
  pheno1 <- pull.pheno(cross1)
  pheno2 <- pull.pheno(cross2)
  correlationData <- vector("list",ncol(pheno1))
  for(x in 1:ncol(pheno1)){
    s <- proc.time()
    r <- apply(pheno2,2,
      function(x,l){
        FUN(x,l,use="pairwise.complete.obs")
      },pheno1[,x])
    
    names(r) <- colnames(pheno2)
    sel <-order(abs(r),decreasing = TRUE)[1:10]
    cat(names(r[sel]),"\n")
    correlationData[[x]] <- r[sel]
    e <- proc.time()
    if(verbose) cat("Item",colnames(pheno1)[x],"max:",max(r,na.rm=T),"min:",min(r,na.rm=T),"took:",round((e[3]-s[3]),digits=2),"Seconds\n",sep=" ")
  }
  correlationData
}
