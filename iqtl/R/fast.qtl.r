# R script for fast QTL via aov function
# (functions are part of the iQTL package)
# Copyright (c) 2010-2012 Danny Arends and Ritsert C. Jansen
# First written Jan, 2011, Last modified Mar, 2012

# Map all phenotypes at a single marker using the aov function
map.fast <- function(marker, phenotypes, conditions, verbose = FALSE){
  st  <- proc.time()
  res <- NULL
  if(is.null(conditions)){
    models  <- aov(as.matrix(phenotypes) ~ marker)
    modelinfo <- summary(models)
    res$qtl <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
    res$eff <- unlist(models$coefficients[2,])
  }else{
    models  <- aov(as.matrix(phenotypes) ~ conditions + marker + conditions:marker)
    modelinfo <- summary(models)
    res$env <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
    res$qtl <- unlist(lapply(modelinfo,"[",2,5),use.names=T)
    res$int <- unlist(lapply(modelinfo,"[",3,5),use.names=T)
    res$eff <- unlist(models$coefficients[2,])
  }
  if(verbose) cat("Done in:", (proc.time()-st)[3], "secs\n")
  res
}

# Create matrices from the per marker modeling results
create.matrix <- function(models, what="qtl", phenonames, genonames, do.log=TRUE){
  res <- matrix(unlist((lapply(models,"[",what))),length(phenonames),length(genonames))
  rownames(res) <- phenonames
  colnames(res) <- genonames
  if(do.log) res <- -log10(res)
  res
}

# Environmental mapping of QTL
# Returns: $qtl, $eff, $env, $int
map.qtl <- function(genotypes, phenotypes, conditions=NULL, n.core=2, verbose = TRUE){
  st  <- proc.time()
  if("parallel" %in% rownames(installed.packages())){
    require("parallel")
    cl <- makeCluster(n.core)
    models <- parApply(cl, genotypes, 2, "map.fast", phenotypes, conditions)
    stopCluster(cl)
  }else{
    models <- apply(genotypes, 2, "map.fast", phenotypes, conditions)  
  }
  res <- NULL
  res$qtl <- create.matrix(models, "qtl", colnames(phenotypes), colnames(genotypes))
  res$eff <- create.matrix(models, "eff", colnames(phenotypes), colnames(genotypes), FALSE)
  if(!is.null(conditions)){
    res$env <- create.matrix(models, "env", colnames(phenotypes), colnames(genotypes))
    res$int <- create.matrix(models, "int", colnames(phenotypes), colnames(genotypes))
  }
  if(verbose)cat("Done in:",(proc.time()-st)[3],"seconds\n")
  res
}

# Main effects permutation vector of maximum QTL
permute.qtl <- function(phenotypes, genotypes, nperm, conditions=NULL, n.core=2, verbose = TRUE){
  indices <- 1:nrow(genotypes)
  permutations <- numeric(nperm)
  for(i in 1:nperm){
    permutations[i] <- max(map.qtl(genotypes[sample(indices),], phenotypes, conditions, n.core)$qtl)
    if(verbose) cat("Permutation", i, "done\n")
  }
  permutations
}
