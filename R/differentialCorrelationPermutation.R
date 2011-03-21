#
#
# differentialCorrelationPermutation.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified feb, 2011
# first written nov, 2010
# 
# R functions to do differential correlation mapping
# Example data here is from C. Elegans and available at request ( Danny.Arends@gmail.com )
# Main bottleneck system memory: (per DifCor object: 3* # of phenotype * #of phenotype)
# Minor bottleneck system cpu: Depending on method, # of phenotype and # of markers
#
# Note: The differentialCorrelation function saves and overwrites the difCor object in the 'output' directory
#

diffCorPermutation <- function(cross, n.perm=10, directory="permutations", verbose=TRUE, ...){
  if(!file.exists(directory)) dir.create(directory)
  for(x in 1:n.perm){
    sl <- proc.time()
    if(verbose) cat("- Starting permutation",x,"/",n.perm,"\n")
		cross$pheno <- cross$pheno[sample(nind(cross)),]
    write.table(diffCorAnalysis(cross, ..., doplot=FALSE, writefile=FALSE, verbose=TRUE),paste(directory,"/Permutation_",x,".txt",sep=""))
    el <- proc.time()
    if(verbose) cat("- Permutation",x,"took:",as.numeric(el[3]-sl[3]),"Seconds.\n")
  }
  difCntPerm <- read.diffCorPermutation(directory)
  invisible(difCntPerm)
}

read.diffCorPermutation <- function(directory="permutations"){
  files <- dir(directory)
  n.perm <- length(files)
  if(n.perm < 1) stop(paste("No permutation files found in:",directory))
  difCntPerm <- vector("list", n.perm)
  for(x in 1:n.perm){
    cat("Trying to read:",paste(directory,"/Permutation_",x,".txt\n",sep=""))
    difCntPerm[[x]]  <- read.table(paste(directory,"/Permutation_",x,".txt",sep=""))
  }
  invisible(difCntPerm)
}

significance.eQCL <- function(difCntPerm){
  maximums <- lapply(difCntPerm,function(x){apply(x,2,max)})
  sorted <- sort(unlist(maximums))
  l <- length(sorted)
  values <- NULL
  valnames <- NULL
  for(x in c(.95,.99,.999)){
    values <- c(values,sorted[l*x])
    valnames <- c(valnames,paste((1-x)*100,"%"))
    cat((1-x)*100,"%\t",sorted[l*x],"\n")
  }
  names(values) <- valnames
  invisible(values)
}

significance.eQCLHotSpot <- function(difCntPerm,significant=212){
  maximums <- lapply(perms,function(x){max(apply(x,1,function(x){sum(x > significant)}))})
  sorted <- sort(unlist(maximums))
  l <- length(sorted)
  values <- NULL
  valnames <- NULL
  for(x in c(.95,.99,.999)){
    values <- c(values,sorted[l*x])
    valnames <- c(valnames,paste((1-x)*100,"%"))
    cat((1-x)*100,"%\t",sorted[l*x],"\n")
  }
  names(values) <- valnames
  invisible(values)
}

