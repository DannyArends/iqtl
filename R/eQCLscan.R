#
# eQCLscan.R
#
# copyright (c) 2010 GBIC : Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Okt, 2011
# first written nov, 2010
# 
# R functions to do eQCL mapping
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

# Counts the number of occurences above threshold in a vector
countVectorThreshold <- function(vector,threshold = 0.5){
  length(which(abs(vector) >= threshold))
}

# Counts the number of occurences above threshold in the eQCL matrix
counteQCLThreshold <- function(eQCL,threshold = 0.5){
  apply(eQCL[[1]],1,countVectorThreshold,threshold)
}

#Print a summary of an eQCL matrix object
summary.eQCL <- function(eQCLMatrix){
  for(s in c(seq(5,50,5),100,200,500,1000)){
    below <- apply(eQCLtMatrix,2,function(x){sum(x>s)==0})
    above <- apply(eQCLMatrix,2,function(x){sum(x>s)!=0})
    cat(s," ",sum(above)," ",sum(below),"\n")
  }
}

#Internal function to correct for differences in mean effect
marker.correct <- function(cross, marker, verbose=FALSE){
  genotypes <- pull.geno(cross)
  phenotypes <- pull.pheno(cross)
  phenotypes <- matrix(as.numeric(unlist(phenotypes)),nrow(phenotypes),ncol(phenotypes))
  m1 <- which(genotypes[,marker]==1)
  m2 <- which(genotypes[,marker]==2)
  oamean <- apply(phenotypes,2,mean,na.rm=T)
  m1mean <- apply(phenotypes[m1,],2,mean,na.rm=T)
  m2mean <- apply(phenotypes[m2,],2,mean,na.rm=T)
  result <- matrix(0,nrow(phenotypes),ncol(phenotypes))
  for(x in m1){ result[x,] <- (oamean-m1mean) }
  for(x in m2){ result[x,] <- (oamean-m2mean) }
  cross$pheno <- as.data.frame(phenotypes + result)
  colnames(cross$pheno) <- phenames(cross)
  cross
}

#Typical analysis routine
#- Load data and remove the NA values from the matrix
#- Detected batch effect and add estimated effect to correct (Could be skipped)
#- Remove QTL effects from the data (Could be skipped)
#Using the Batch and QTL correction matrix
yeast.analysis.eQTL <- function(directory="E:/GBIC/Bruno/Differential Correlations/Yeast2"){
  require(qtl)
  require(iqtl)
  memory.limit(4000)
  setwd(directory)
  bremcross <- read.cross("csvr",file="yeast_brem_cross.csv",geno=c(0,1),colClasses="numeric")
  bremcross <- convert2riself(bremcross)
  newpheno <- as.data.frame(as.matrix(apply(bremcross$pheno,2,as.numeric),nrow(bremcross$pheno),ncol(bremcross$pheno)))
  bremcross$pheno <- newpheno
  #phenotypes <- pull.pheno(bremcross)
  #phenotypevariance <- apply(bremcross$pheno,2,var)
  #genes <- which(phenotypevariance >= 0)
  #bremcross$pheno <- phenotypes[,genes]
  #batchcorrection <- read.table("batchmatrix.txt",header=T)
  #bremcross$pheno <- bremcross$pheno-batchcorrection
  #qtlcorrection <- read.table("qtlmatrix.txt",header=T)
  #bremcross$pheno <- bremcross$pheno+qtlcorrection
  MyeQCLMatrix <- eQCLscan(bremcross,0,0.4,50,pheno.col=43:282,directory="WithQTL",doplot=T,writefile=T)
  MyeQCLCnts <- eQCLpermute(bremcross, eQCLThreshold=0.4)
}

#Test (and time) the eQCLscan routine
#Using the Batch and QTL correction matrix
test.eQTLperformance <- function(directory="e:/gbic/bruno/differential correlations/yeast2"){
  require(qtl)
  setwd(directory)
  cross <- read.cross("csvr",file="yeast_brem_cross.csv",geno=c(0,1))
  cross <- convert2riself(bremcross)
  batchcorrection <- read.table("batchmatrix.txt",header=T)
  qtlcorrection <- read.table("qtlmatrix.txt",header=T)
  cross$pheno <- cross$pheno-batchcorrection
  cross$pheno <- cross$pheno+qtlcorrection
  times <- NULL
  for(x in seq(100,1000,50)){
    s <- proc.time()
    eQCLscan(cross,1:x)
    e <- proc.time()
    times <- c(times,as.numeric(e[3]-s[3]))
  }
  names(times) <- paste("items_",seq(100,1000,50),sep="")
  times
}

#Create a correlation matrix (so we can more easily use a lineair vectory when using apply)
correlationmatrix <- function(x,expressions,method="pearson"){
  cor(expressions[x,], method = method, use="pairwise.complete.obs")
}

#eQCL.scanmarker routine - Optimized using snow and 2 cores
#Does a single marker return a list with:
# [[1]] eQCLmatrix 
# [[2]] Correlation matrix individuals with genotype 1
# [[3]] Correlation matrix individuals with genotype 2
# [[4]] A vector of counted eQCL scores for each phenotype above the eQCLThreshold
#Note: Also saves the object to: output/eQCL<marker>.Rdata
eQCL.scanmarker <- function(cross, marker, eQCLThreshold=0.25, method="pearson", directory="output", saveRdata=FALSE,snow=TRUE){
  if(is.missing(cross)) stop("No cross object")
  if(is.missing(marker)) warning("Doing all markers, be warned that this will take up A LOT of RAM memory")
  require(snow)
  expressions <- matrix(unlist(pull.pheno(cross)),dim(pull.pheno(cross))[1],dim(pull.pheno(cross))[2])
  colnames(expressions) <- colnames(pull.pheno(cross))
  genotypes <- pull.geno(cross)
  markerName <- markernames(cross)[marker]
  
  work <- vector("list",2)
  work[[1]] <- which(genotypes[,marker]==1)
  work[[2]] <- which(genotypes[,marker]==2)
  
  if(snow){
    cpu_cluster <- makeCluster(c("localhost","localhost"))
    results <- parLapply(cpu_cluster,work,correlationmatrix,expressions=expressions,method=method)
    stopCluster(cpu_cluster)
  }else{
    results <- lapply(work,correlationmatrix,expressions=expressions,method=method)
  }
  
  eQCL <- vector("list",4)
  eQCL[[1]] <- 0.5*((sign(results[[1]]) * results[[1]]^2) -  (sign(results[[2]]) * results[[2]]^2))
  eQCL[[2]] <- results[[1]]
  eQCL[[3]] <- results[[2]]
  eQCL[[4]] <- counteQCLThreshold(eQCL,eQCLThreshold)
  
  traitnames <- colnames(expressions)
  colnames(eQCL[[1]]) <- traitnames
  rownames(eQCL[[1]]) <- traitnames
  colnames(eQCL[[2]]) <- traitnames
  rownames(eQCL[[2]]) <- traitnames
  colnames(eQCL[[3]]) <- traitnames
  rownames(eQCL[[3]]) <- traitnames
  
  attr(eQCL,"marker") <- markerName
  attr(eQCL,"phenotypes") <- traitnames
  
  if(saveRdata) save(eQCL,file=paste(directory,"/eQCL",marker,".Rdata",sep=""))
  
  eQCL
}

#Extract a list of traits which show eQCL at a marker
#Returns a list of traits which show eQCL
eQCL.targets <- function(eQCL, eQCLThreshold=0.25, p.value = 0, verbose=TRUE){
  if(is.missing(eQCL)) stop("No eQCL object")
  aboveThreshold <- counteQCLThreshold(eQCL, eQCLThreshold)
  eQCLrelated <- which(aboveThreshold > significant)
  cat("Traits showing eQCL: ",length(names(eQCLrelated)),"\n")
  if(length(names(eQCLrelated)) > 1){
    eQCLMatrix <- eQCL[[1]][names(eQCLrelated),]
  }else{
    eQCLMatrix <- t(as.matrix(eQCL[[1]][names(eQCLrelated),],ncol(eQCL[[1]]),1))
    rownames(eQCLMatrix) <- names(eQCLrelated)
  }
  out <- vector("list",length(names(eQCLrelated)))
  cnt <- 1
  for(x in names(eQCLrelated)){
    out[[cnt]] <- colnames(eQCLMatrix)[which(abs(eQCLMatrix[x,]) > eQCLThreshold)]
    attr(out[[cnt]],"name") <- x
    if(verbose) cat(x," ",length(out[[cnt]]),"\n")
    cnt <- cnt+1
  }
  out
}

#Scale down phenotypes
#We throw away phenotypes which have less variation then the minimumvariance
scaledownPhenotypes <- function(cross,minimumvariance = 0.1,verbose=FALSE){
  phenotypes <- pull.pheno(cross)
  phenotypevariance <- apply(cross$pheno,2,var,na.rm=T)
  genes <- which(phenotypevariance >= minimumvariance)
  cross$pheno <- phenotypes[,genes]
  if(verbose) cat(nphe(cross),"/",ncol(phenotypes)," remaining\n")
  cross
}

#Get the names of the probes that are significant
eQCL.names <- function(eQCLMatrix, significance=212){
  names(which(apply(eQCLMatrix,2,function(x){max(x)>significance})))
}

annotate.group <- function(group,annotation){
  annotation[which(annotation[,1] %in% group),2]
}

#Main routine to do the entire analysis
#Note: Does all the markers one by one (optimized to use 2 cores)
#Note: The eQCL object in memory is very large
#Note: Based on the amount of traits and markers this could take a LONG time
eQCLscan <- function(cross, marker.range, eQCLThreshold=0.25, significant = 5, method="pearson", directory="output", doplot=FALSE, writefile=FALSE, saveRdata=FALSE, snow=TRUE, verbose=TRUE){
  if(is.missing(cross)) stop("No cross object")
  s <- proc.time()
  if(doplot && !file.exists(directory)) dir.create(directory)
  if(verbose) cat("Analysis of ",ncol(cross$pheno)," traits at ",sum(nmar(cross))," markers\n")
  if(is.missing(marker.range)){
    warning("Selecting all markers this will take a LONG time")
    marker.range <- 1:sum(nmar(cross))
  }
  eQCLMatrix <- NULL
  
  for(marker in marker.range){
    sl <- proc.time()
    if(snow){
      p_usage <- gc()[2,3]
      n_usage <- gc()[2,3]
      while(n_usage < p_usage){
        p_usage = n_usage
        n_usage <- gc()[2,3]
        if(verbose) cat("GCloop ",n_usage," ",p_usage,"\n")
      }
    }
    results <- eQCL.scanmarker(cross, marker, eQCLThreshold, method,directory,saveRdata,snow)
    if(doplot){
      png(paste(directory,"/",marker,"_",markernames(cross)[marker],".jpg",sep=""))
      plot.eQCL(results, eQCLThreshold, significant)
      dev.off()
    }
    
    eQCLMatrix <- rbind(eQCLMatrix,results[[4]])
    
    results <- NULL
    el <- proc.time()
    if(verbose){
      cat("Marker ",marker,"/ [",min(totmarkers),"..",max(totmarkers),"]\n")
      cat("Time: ",as.numeric(el[3]-sl[3]),", so far:",as.numeric(el[3]-s[3]),"seconds.\n")
    }
  }
  cat("Analysis took: ",as.numeric(el[3]-s[3]),"seconds\n")
  rownames(eQCLMatrix) <- markernames(cross)[totmarkers]
  colnames(eQCLMatrix) <- phenames(cross)
  if(writefile) write.table(eQCLMatrix,file="eQCLMatrix.txt",sep="\t")
  eQCLMatrix
}

