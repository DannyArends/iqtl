#
# eQCLpermute.R
#
# copyright (c) 2010 GBIC : Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Okt, 2011
# first written nov, 2010
# 
# R functions to do permutation on eQCL mapping, in order to determin significance
# Example datasets are available at request ( Danny.Arends@gmail.com )
#  C. Elegans
#  Yeast
#  Brassics Napus/Rapa
#  Arabidopsis
#  
# Main bottleneck system memory: (per eQCL object: 3* # of phenotype * # of phenotype)
# Minor bottleneck system cpu: Depending on method, # of phenotype and # of markers
#
# Note: The eQCL function saves and overwrites the eQCL object in the 'output' directory
#

eQCLpermute <- function(cross, n.perm=10, directory="permutations", verbose=TRUE, ...){
  if(!file.exists(directory)) dir.create(directory)
  for(x in 1:n.perm){
    sl <- proc.time()
    if(verbose) cat("- Starting permutation",x,"/",n.perm,"\n")
    cross$pheno <- cross$pheno[sample(nind(cross)),]
    write.table(eQCLscan(cross, ..., doplot=FALSE, writefile=FALSE, verbose=TRUE),paste(directory,"/perm_",x,".txt",sep=""))
    el <- proc.time()
    if(verbose) cat("- Permutation",x,"took:",as.numeric(el[3]-sl[3]),"seconds.\n")
  }
  eQCLperm <- load.eQCLpermute(directory)
  invisible(eQCLperm)
}

#TODO: Use a regexp to select the perm_ files
load.eQCLpermute <- function(directory="permutations"){
  files <- dir(directory)
  n.perm <- length(files)
  if(n.perm < 1) stop(paste("No permutation files found in:",directory))
  eQCLperm <- vector("list", n.perm)
  for(x in 1:n.perm){
    cat("Trying to read:",paste(directory,"/perm_",x,".txt\n",sep=""))
    eQCLperm[[x]]  <- read.table(paste(directory,"/perm_",x,".txt",sep=""))
  }
  invisible(eQCLperm)
}

significance.eQCL <- function(eQCLperm, alpha, verbose = FALSE){
  if(is.missing(alpha)){
    maximums <- lapply(eQCLperm,function(x){apply(x,2,max)})
  }else{
    maximums <- lapply(eQCLperm,function(x){max(apply(x,1,function(x){sum(x > alpha)}))})
  }
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

