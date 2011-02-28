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

#Typical analysis routine
#- Load data and remove the NA values from the matrix
#- Detected batch effect and add estimated effect to correct (Could be skipped)
#- Remove QTL effects from the data (Could be skipped)
#Using the Batch and QTL correction matrix
analysis.differentialCorrelation <- function(){
	require(qtl)
  require(iqtl)
  memory.limit(3000)
	setwd("e:/gbic/bruno/differential correlations/yeast2")
	bremcross <- read.cross("csvr",file="yeast_brem_cross.csv",geno=c(0,1))
  bremcross <- convert2riself(bremcross)
  phenotypes <- pull.pheno(bremcross)
  phenotypevariance <- apply(bremcross$pheno,2,var)
  genes <- which(phenotypevariance >= 0)
  bremcross$pheno <- phenotypes[,genes]
  batchcorrection <- read.table("batchmatrix.txt",header=T)
  bremcross$pheno <- bremcross$pheno-batchcorrection
  qtlcorrection <- read.table("qtlmatrix.txt",header=T)
  bremcross$pheno <- bremcross$pheno+qtlcorrection
  MydifCntMatrix <- diffCorAnalysis(bremcross,0.01,0.4,10,directory="WithOutQTL")
  
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
#
#Note: Also saves the object to: output/difCor<marker>.Rdata
differentialCorrelation <- function(cross, marker, difCorThreshold=0.5, method="pearson", directory="output", saveRdata=FALSE){
  require(snow)
  expressions <- matrix(unlist(pull.pheno(cross)),dim(pull.pheno(cross))[1],dim(pull.pheno(cross))[2])
  colnames(expressions) <- colnames(pull.pheno(cross))
  genotypes <- pull.geno(cross)
  markerName <- markernames(cross)[marker]
  
  work <- vector("list",2)
  work[[1]] <- which(genotypes[,marker]==1)
  work[[2]] <- which(genotypes[,marker]==2)
  
  cpu_cluster <- makeCluster(c("localhost","localhost"))
  results <- parLapply(cpu_cluster,work,correlationmatrix,expressions=expressions,method=method)
  stopCluster(cpu_cluster)
  
  difCor <- vector("list",4)
  difCor[[1]] <- (sign(results[[1]]) * results[[1]]^2) -  (sign(results[[2]]) * results[[2]]^2)
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

#Heatmap the output of a difCor object
#Returns the hclust object used to order the traits shown in the heatmap
plotDifCor <- function(difCor, difCorThreshold=0.5, significant = 0, ...){
  aboveThreshold <- countDifCorThreshold(difCor, difCorThreshold)
  difCorrelated <- which(aboveThreshold > significant)
  ccorclass1 <- difCor[[2]][difCorrelated,difCorrelated]
  ccorclass2 <- difCor[[3]][difCorrelated,difCorrelated]
  if(nrow(ccorclass1) >= 2){
    ordering <- hclust(dist(ccorclass1))$order
    ccorclass1 <- ccorclass1[ordering,ordering]
    ccorclass2 <- ccorclass2[ordering,ordering]
    upper <- t(upper.tri(ccorclass1))*ccorclass1
    lower <- t(lower.tri(ccorclass2))*ccorclass2
    op <- par(mar=c(8,8,4,2)+0.1)
    genotypes <- upper.tri(matrix(0,nrow(ccorclass1),nrow(ccorclass1)))
    for(x in 1:nrow(genotypes)){
      genotypes[x,x] <- NA
    }
    image(seq(0,nrow(ccorclass1)),seq(0,nrow(ccorclass1)),genotypes,ylab="",xlab="",yaxt="n",xaxt="n",col=c(rgb(255,0,0,50,maxColorValue=255),rgb(255,255,0,50,maxColorValue=255)),main=paste("DifCor at marker",attr(difCor,"marker")),...)
    image(seq(0,nrow(ccorclass1)),seq(0,nrow(ccorclass1)),lower+upper,breaks=c(-1,-0.75,-0.5,0.5,0.75,1),col=c("blue","lightblue",rgb(255,255,255,155,maxColorValue=255),"lightgreen","green"),add=T,...)
    axis(2,at=seq(0.5,nrow(ccorclass1)),labels=colnames(ccorclass1),las=1)
    axis(1,at=seq(0.5,nrow(ccorclass1)),labels=colnames(ccorclass1),las=3)
    if(nrow(genotypes) < 50){
      grid(nrow(ccorclass1),nrow(ccorclass1),lwd=1,lty=1,col="black")
    }else{
      cat("INFO: ",attr(difCor,"marker"),"- 50+ traits, plot grid disabled, Add with: grid(",nrow(genotypes),",",nrow(genotypes),",lwd=1,lty=1,col=\"black\")\n")
    }
    abline(v=0)
    abline(h=nrow(ccorclass1))
    ordering
  }else{
    warning("No phenotype shows differential correlation with more than: ",significant," other phenotypes at difCorThreshold: ",difCorThreshold,"\n")
  }
}

#Returns the matrix with differential correlations counts.
#Counts are: differential correlations above a threshold
# rows: Phenotypes
# cols: Markers
#difCntMatrix <- function(difCor, directory="output", difCorThreshold=0.5){
#  cntMatrix <- NULL
#  markerNames <- NULL
#  for(marker in 1:length(dir(directory))){
#    load(paste(directory,"/difCor",marker,".Rdata",sep=""))
#    cat(attr(difCor,"marker"),"\n")
#    cntMatrix <- rbind(cntMatrix,countDifCorThreshold(difCor,difCorThreshold))
#    markerNames <- c(markerNames,attr(difCor,"marker"))
#    phenotypeNames <- attr(difCor,"phenotypes")
#  }
#  rownames(cntMatrix) <- markerNames
#  colnames(cntMatrix) <- phenotypeNames
#  cntMatrix
#}

#Heatmap the output of a difCnt object
imageDifCnt <- function(difCnt, cluster=FALSE){

}

#Plot a Differential Correlation Count profile of a single phenotype
#Optionally add a scanone object to overlay the QTL profile
plotDifCntProfile <- function(cross, difCntMatrix, pheno.col=1, addQTL=FALSE){
  difCorCntProfile <- lodscorestoscanone(cross,difCntMatrix[,pheno.col],traitnames = "difCorCnt")
  if(addQTL) qtlscan <- scanone(cross, pheno.col=pheno.col)
  dtmax <- max(difCorCntProfile[,3])
  if(addQTL) qtlmax <- max(qtlscan[,3])
  if(addQTL) difCorCntProfile[,3] <- difCorCntProfile[,3]*(qtlmax/dtmax)
  if(addQTL){
    plot(qtlscan,difCorCntProfile,y=c(0,1.7*qtlmax),col=c("red","green"),main=pheno.col)
    legend("topright",c("scanone","difCorCount"),lwd=1,col=c("red","green"))
    axis(4,at=seq(0,1.7*qtlmax,1),round((dtmax/qtlmax) * seq(0,1.7*qtlmax,1),1))
  }else{
    plot(difCorCntProfile,y=c(0,1.7*dtmax),col="black",main=pheno.col,ylab="DifCorCnt")
    legend("topright",c("difCorCount"),lwd=1,col=c("black"))
  }
  if(addQTL) difCorCntProfile[,3] <- difCorCntProfile[,3]*(dtmax/qtlmax)
  difCorCntProfile
}

plotDifCorAtMarker <- function(cross, difCntMatrix, marker="YBR008C_211", significant=5, lodthreshold=5){
  signdifcor <- names(which(difCntMatrix[marker,] > significant))
  cat("- Starting QTL scan of",length(signdifcor),"\n")
  res <- scanone(cross,pheno.col=signdifcor)
  cat("- Start of",length(which(apply(res[,3:ncol(res)],2,max) > lodthreshold)),"plots\n")
  if(length(which(apply(res[,3:ncol(res)],2,max) > lodthreshold)) > 0){
    for(x in names(which(apply(res[,3:ncol(res)],2,max) > lodthreshold))){
      plotDifCntProfile(cross,difCntMatrix,x)
    }
  }
  invisible(res)
}

#Plot the change in Correlation profile of a single phenotype
#Not a priority
plotDifCorNetwork <- function(difCor, difCorThreshold=0.5, significant = 0){

}


#Extract a list of traits which show differential correlation
#Returns a list of traits which show differential correlation with more then the user specified significant # of genes
difCorTargets <- function(difCor, difCorThreshold=0.5, significant = 0, verbose=TRUE){
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
    attr(out,"name") <- x
    if(verbose) cat(x," ",length(out[[cnt]]),"\n")
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

#Main routine to do the entire analysis
#Get rid of traits with no/low expression variation
#Note: Does all the markers one by one (could be optimized to use 2 cores, however the difCor object in memory is large
#Note: Based on the amount of traits and markers this could take a LONG time
diffCorAnalysis <- function(cross, minimumvariance=0.03, difCorThreshold=0.4, significant = 5, method="pearson", directory="output", verbose=TRUE){
  s <- proc.time()
  if(!file.exists(directory)) dir.create(directory)
  cross <- scaledownPhenotypes(cross, minimumvariance,verbose)
  if(verbose) cat("Analysis of ",ncol(cross$pheno)," traits at ",sum(nmar(cross))," markers\n")
  
  totmarkers <- sum(nmar(cross))
  difCountMatrix <- NULL
  
  for(marker in 1:totmarkers){
    sl <- proc.time()
    results <- differentialCorrelation(cross, marker, difCorThreshold, method,directory)
    png(paste(directory,"/",marker,"_",markernames(cross)[marker],".jpg",sep=""))
      plotDifCor(results, difCorThreshold, significant)
    dev.off()
    
    difCountMatrix <- rbind(results[[4]],difCountMatrix)
    
    results <- NULL
    gc()
    el <- proc.time()
    if(verbose){
      cat("Marker ",marker,"/",totmarkers,"took: ",as.numeric(el[3]-sl[3]),"Seconds.\n")
    }
  }
  cat("Analysis took: ",as.numeric(el[3]-s[3]),"Seconds\n")
  rownames(difCountMatrix) <- markernames(cross)
  colnames(difCountMatrix) <- phenames(cross)
  write.table(difCountMatrix,file="difCountMatrix.txt",sep="")
  difCountMatrix
}
