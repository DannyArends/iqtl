#####################################################################
 # 
 # scanRF.R
 #
 # Copyright (c) 2009 Danny Arends
 # Last modified Dec, 2009
 # First written Mrt, 2009
 #
 # scanRF
 # Contains: scanRF - 	N-Dimensional phenotype additive interactions scanning methode for QTL mapping using the randomForest algorithm
 #			 			          made as a training excersice in QTL mapping in RIL/BC and F2 crosses
 #						          Usefull for fast lowlevel (markerwise) data exploration for QTL mapping experiments
 # Depends: R/QTL package
 # 			randomForest
 # 
 #
#####################################################################

#unittest.scanRF <- function(){
#	require(qtl)
#	require(randomForest)
#	data(listeria)
#	data(hyper)	
#	aaa <- scanRF(listeria)					#Listeria dataset Phenotype1
#	aaa	<- scanRF(hyper)					#Hypertensia dataset Phenotype1
#}

mapRF <- function(data,model){
	require(randomForest)
	temp <- cbind(data$phenotype,t(data$genotypes))
	colnames(temp)[1] <- "phenotype"
	rf <- randomForest(phenotype ~ ., data=temp,importance=T,proximity=T)
	ret <- -10*log10(1-(abs(importance(rf)[,1])/abs(max(importance(rf)[,1])+1)))
  ret
}



scanRF <- function(cross=NULL, pheno.col=1,plot=FALSE , compare=FALSE, ...){
	require(qtl)
	require(randomForest)
	cross <- qtl::fill.geno(cross)

	#Handle phenotype from R/QTL
	pheno <- NULL
	
	if (length(pheno.col) > 1){
		pheno <- cross$pheno[,pheno.col]
		n.traits <- ncol(pheno)
	}else{
		if(qtl::nphe(cross) > pheno.col && pheno.col >0){
			pheno <- qtl::pull.pheno(cross)[pheno.col]
			n.traits <- 1
		}else{
			stop("ERROR: Wrong phenotype..") 
		}
	}
	#Handle genotypes
	data <- NULL
  chr <- NULL
  map <- NULL
	for(i in 1:length(cross$geno)){
		#cat("Chromosome: ",i,"\n")
		data <- cbind(data,cross$geno[[i]]$data)
    chr <- c(chr,rep(i,length(cross$geno[[i]]$map)))
    map <- c(map,as.double(cross$geno[[i]]$map))
	}
	names <- NULL
	#Name the phenotypes
	xnam <- NULL
	for(i in 1:n.traits){
		xnam <- c(xnam,paste("P", pheno.col[i], sep=""))
	}
	names <- xnam
	for(i in 1:dim(data)[2]){
		#Rename markers so we don't have naming problems
		names <- c(names,paste("M",i,sep=""))
	}
	#cat(names,"\n")
	data <- cbind(pheno,data)
  oldnames <- colnames(data)
  oldnames <- oldnames[-1]
	colnames(data) <- names
	#remove all missing phenotypes (we filled in the geno's)
	data <- na.exclude(data)
	if(n.traits==1){
		howmany <-1
	}else{
		howmany <- n.traits + 1
	}
	xnames <- NULL
	for(i in 1:(howmany)){
		if(i == howmany){
			xnames = xnam
			temp <- data
		}else{
			xnames <- c(paste("P", pheno.col[i], sep=""))
			temp <- data[,-which(!(xnam==xnames))]
		}
		rf <- randomForest((fmla <- as.formula(paste(paste(xnames, collapse= "+"), " ~ ."))) , data=temp, importance=T,proximity=T, ...)
	
		resrf <- -10*log10(1-(abs(importance(rf)[,1])/abs(max(importance(rf)[,1])+1)))
		resrf <- cbind(resrf,-10*log10(1-(abs(importance(rf)[,2])/abs(max(importance(rf)[,2])+1))))
	}
	resrf <- cbind(chr,map,resrf)
  rownames(resrf) <- oldnames
  colnames(resrf) <- c("chr","pos","lod","lod2")
  resrf <- as.data.frame(resrf)
  class(resrf) <- c("scanone",class(resrf))
  if(plot){
    if(compare){
      res <- scanone(cross,pheno.col=pheno.col)
      plot(resrf,resrf[,-3],res,main="QTLscan using randomforest",sub="comparing back to scanone",col=c("Black","Gray","Blue"),lwd=c(2,1,2))
      legend("topright",c("scanRF importance measure 1","scanRF importance measure 2","scanone"),col=c("Black","Gray","Blue"),lwd=c(2,1,2))
    }else{
      plot(resrf,lodcolum=c(1,2),main="QTLscan using randomforest",col=c("Black","Gray"),lwd=c(2,1))
      legend("topright",c("scanRF importance measurement 1","scanRF importance measurement 2"),col=c("Black","Gray"),lwd=c(2,1))
    }
  }
  resrf
}
