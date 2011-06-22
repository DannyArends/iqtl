#
#
# differentialCorrelation.R
#
# copyright (c) 2010 Bruno Tesson and Danny Arends
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

# Counts the number of occurences above threshold in a vector
countVectorThreshold <- function(vector,threshold = 0.5){
  length(which(abs(vector) >= threshold))
}

# Counts the number of occurences above threshold in the difCor matrix
countDifCorThreshold <- function(difCor,threshold = 0.5){
  apply(difCor[[1]],1,countVectorThreshold,threshold)
}

summaryDifCnt <- function(difCntMatrix){
  for(s in c(seq(5,50,5),100,200,500,1000)){
    not_significant <- apply(difCntMatrix,2,function(x){sum(x>s)==0})
    significant <- apply(difCntMatrix,2,function(x){sum(x>s)!=0})
    cat(s," ",sum(not_significant)," ",sum(significant),"\n")
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
analysis.differentialCorrelation <- function(){
	require(qtl)
  require(iqtl)
  memory.limit(4000)
	setwd("E:/GBIC/Bruno/Differential Correlations/Yeast2")
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
  MydifCntMatrix <- diffCorAnalysis(bremcross,0,0.4,50,pheno.col=43:282,directory="WithQTL",doplot=T,writefile=T)
  MydifPermCnts <- diffCorPermutation(bremcross, minimumvariance=0, difCorThreshold=0.4)
  #Reload it form disk
  #setwd("e:/gbic/bruno/differential correlations/yeast2")
  #MydifCntMatrix <- difCntMatrix()
}

#Test (and time) the differentialCorrelation routine
#Using the Batch and QTL correction matrix
test.differentialCorrelation <- function(){
	require(qtl)
	setwd("e:/gbic/bruno/differential correlations/yeast2")
	bremcross <- read.cross("csvr",file="yeast_brem_cross.csv",geno=c(0,1))
  bremcross <- convert2riself(bremcross)
  batchcorrection <- read.table("batchmatrix.txt",header=T)
  qtlcorrection <- read.table("qtlmatrix.txt",header=T)
  cross$pheno <- cross$pheno-batchcorrection
  cross$pheno <- cross$pheno+qtlcorrection
	times <- NULL
	for(x in seq(100,1000,50)){
		s <- proc.time()
		differentialCorrelation(bremcross,1:x)
		e <- proc.time()
		times <- c(times,as.numeric(e[3]-s[3]))
	}
	names(times) <- paste("items_",seq(100,1000,50),sep="")
  times
}

#Create a correlation matrix (so we can easily apply
correlationmatrix <- function(x,expressions,method="pearson"){
	cor(expressions[x,], method = method, use="pairwise.complete.obs")
}

#Differential correlation routine - Optimized using snow and 2 cores
#Does a Single marker return a list with:
# [[1]] difCorMatrix 
# [[2]] Correlation matrix individuals with genotype 1
# [[3]] Correlation matrix individuals with genotype 2
# [[4]] A vector of counted scores for each phenotype above the difCorThreshold
#Note: Also saves the object to: output/difCor<marker>.Rdata
differentialCorrelation <- function(cross, marker, difCorThreshold=0.25, method="pearson", directory="output", saveRdata=FALSE,snow=TRUE){
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
  
  difCor <- vector("list",4)
  difCor[[1]] <- 0.5*((sign(results[[1]]) * results[[1]]^2) -  (sign(results[[2]]) * results[[2]]^2))
  difCor[[2]] <- results[[1]]
  difCor[[3]] <- results[[2]]
  difCor[[4]] <- countDifCorThreshold(difCor,difCorThreshold)
  
  traitnames <- colnames(expressions)
  colnames(difCor[[1]]) <- traitnames
  rownames(difCor[[1]]) <- traitnames
  colnames(difCor[[2]]) <- traitnames
  rownames(difCor[[2]]) <- traitnames
  colnames(difCor[[3]]) <- traitnames
  rownames(difCor[[3]]) <- traitnames
  
  attr(difCor,"marker") <- markerName
  attr(difCor,"phenotypes") <- traitnames
  
  if(saveRdata) save(difCor,file=paste(directory,"/difCor",marker,".Rdata",sep=""))
  
  difCor
}

#Extract a list of traits which show differential correlation
#Returns a list of traits which show differential correlation with more then the user specified significant # of genes
difCorTargets <- function(difCor, difCorThreshold=0.25, significant = 0, verbose=TRUE){
  aboveThreshold <- countDifCorThreshold(difCor, difCorThreshold)
  difCorrelated <- which(aboveThreshold > significant)
  cat("Traits showing difCor: ",length(names(difCorrelated)),"\n")
  if(length(names(difCorrelated)) > 1){
    difCorMatrix <- difCor[[1]][names(difCorrelated),]
  }else{
    difCorMatrix <- t(as.matrix(difCor[[1]][names(difCorrelated),],ncol(difCor[[1]]),1))
    rownames(difCorMatrix) <- names(difCorrelated)
  }
  out <- vector("list",length(names(difCorrelated)))
  cnt <- 1
  for(x in names(difCorrelated)){
    out[[cnt]] <- colnames(difCorMatrix)[which(abs(difCorMatrix[x,]) > difCorThreshold)]
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
eQCL.significant <- function(difCntMatrix, significance=212){
  names(which(apply(difCntMatrix,2,function(x){max(x)>significance})))
}

annotate.group <- function(group,annotation){
  annotation[which(annotation[,1] %in% group),2]
}

#Main routine to do the entire analysis
#Get rid of traits with no/low expression variation
#Note: Does all the markers one by one (optimized to use 2 cores)
#Note: The difCor object in memory is very large
#Note: Based on the amount of traits and markers this could take a LONG time
diffCorAnalysis <- function(cross, difCorThreshold=0.25, significant = 5, marker.col=NULL, method="pearson", directory="output", doplot=FALSE, writefile=FALSE, saveRdata=FALSE, snow=TRUE, verbose=TRUE){
  s <- proc.time()
  if(doplot && !file.exists(directory)) dir.create(directory)
  #cross <- scaledownPhenotypes(cross, minimumvariance,verbose)
  if(verbose) cat("Analysis of ",ncol(cross$pheno)," traits at ",sum(nmar(cross))," markers\n")
  if(!is.null(marker.col)){
    totmarkers <- marker.col
  }else{
    totmarkers <- 1:sum(nmar(cross))
  }
  difCountMatrix <- NULL
  
  for(marker in totmarkers){
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
    results <- differentialCorrelation(cross, marker, difCorThreshold, method,directory,saveRdata,snow)
    if(doplot){
      png(paste(directory,"/",marker,"_",markernames(cross)[marker],".jpg",sep=""))
        plotDifCor(results, difCorThreshold, significant)
      dev.off()
    }
    
    difCountMatrix <- rbind(difCountMatrix,results[[4]])
    
    results <- NULL
    el <- proc.time()
    if(verbose){
      cat("Marker ",marker,"/ [",min(totmarkers),"..",max(totmarkers),"] took: ",as.numeric(el[3]-sl[3]),", so far:",as.numeric(el[3]-s[3]),"Seconds.\n")
    }
  }
  cat("Analysis took: ",as.numeric(el[3]-s[3]),"Seconds\n")
  rownames(difCountMatrix) <- markernames(cross)[totmarkers]
  colnames(difCountMatrix) <- phenames(cross)
  if(writefile) write.table(difCountMatrix,file="difCountMatrix.txt",sep="\t")
  difCountMatrix
}
