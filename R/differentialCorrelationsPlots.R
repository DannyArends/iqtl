#
#
# differentialCorrelationPlots.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Jun, 2011
# first written nov, 2010
# 
# Plotting routines for differentialCorrelation Analysis
#

#Heatmap the output of a difCor object
#Returns the hclust object used to order the traits shown in the heatmap
plotDifCor <- function(difCor, difCorThreshold=0.5, significant = 0, ...){
  aboveThreshold <- countDifCorThreshold(difCor, difCorThreshold)
  difCorrelated <- which(aboveThreshold > significant)
  if(length(difCorrelated) <= 1){
     warning("No phenotype shows differential correlation with more than: ",significant," other phenotypes at difCorThreshold: ",difCorThreshold,"\n")
     return()
  }
  ccorclass1 <- difCor[[2]][difCorrelated,difCorrelated]
  ccorclass2 <- difCor[[3]][difCorrelated,difCorrelated]
  if(nrow(ccorclass1) >= 2){
    ordering <- hclust(dist(ccorclass1))
    ccorclass1 <- ccorclass1[ordering$order,ordering$order]
    ccorclass2 <- ccorclass2[ordering$order,ordering$order]
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

#Heatmap the output of a difCor object (almost the same as above)
imageDifCor <- function(difcor, peekheight=5){
  selection <- which(difcor[[4]] > peekheight)
  if(length(selection) > 0){
    g1 <- difcor[[2]][selection,selection]
    g2 <- difcor[[3]][selection,selection]
    ordering <- difcor[[1]][selection,selection]
    clustering <- hclust(dist(ordering))

    upper <- upper.tri(g1)*sign(g1)*(g1[clustering$order,clustering$order]^2)
    lower <- lower.tri(g2)*sign(g2)*(g2[clustering$order,clustering$order]^2)
    colorz <- c("red","white","blue")
    heatmap(t(upper+lower),col=colorz,breaks=c(-1,-0.25,0.25,1),Colv=NA,Rowv=NA,scale="none",main=paste("Correlation at: ",attr(difCor,"marker")))
    clustering
  }else{
    plot(1:10)
    return(NULL)
  }
}

#Heatmap the output of a difCnt object
imageDifCnt <- function(difCnt, significant=100, cluster=FALSE){
  if(cluster){
    heatmap(t(difCnt[,which(apply(difCnt,2,function(x){max(x) > significant}))]),Colv=NA,breaks=c(0,100,10000),scale="none",col=c("white","black"))
  }else{
    heatmap(t(difCnt[,which(apply(difCnt,2,function(x){max(x) > significant}))]),Colv=NA, Rowv=NA,breaks=c(0,100,10000),scale="none",col=c("white","black"))
  }
}

#plotDifCntProfile Plot a Differential Correlation Count profile of a single phenotype
# Optionally add a scanone object to overlay the QTL profile
plotDifCntProfile <- function(cross, difCntMatrix, pheno.col=1,significanceThresholds=NULL, addQTL=FALSE, ...){
  difCorCntProfile <- lodscorestoscanone(cross,difCntMatrix[,pheno.col],traitnames = "difCorCnt")
  if(addQTL){
    cross <- calc.genoprob(cross)
    qtlscan <- scanone(cross, pheno.col=pheno.col)
  }
  dtmax <- max(difCorCntProfile[,3])
  if(addQTL) qtlmax <- max(qtlscan[,3],10)
  if(addQTL) difCorCntProfile[,3] <- difCorCntProfile[,3]*(qtlmax/dtmax)
  if(addQTL){
    op <- par(mar=c(5, 4, 4, 5) + 0.1)
    plot(qtlscan,difCorCntProfile,y=c(0,1.7*qtlmax),col=c("red","black"),lty=c(1,1),lwd=c(3,2),main=pheno.col,...)
    axis(4,at=seq(0,1.7*qtlmax,1),round((dtmax/qtlmax) * seq(0,1.7*qtlmax,1),1))
    legend("topright",c("scanone","difCorCount"),lty=c(1,1),lwd=c(3,2),col=c("red","black"))
  }else{
    plot(difCorCntProfile,y=c(0,1.7*dtmax),col="black",main=pheno.col,ylab="DifCorCnt",...)
    if(!is.null(significanceThresholds)){
      colorz <- c("red","orange","green")
      i <- 1
      for(x in significanceThresholds){
        abline(h=x,col=colorz[i],lty=3,lwd=2)
        i <- i+1
      }
      legend("topright",c("difCorCount",names(significanceThresholds)),lwd=2,col=c("black",colorz),lty=c(1,3,3,3))
    }else{
      legend("topright",c("difCorCount"),lwd=1,col=c("black"))
    }
  }
  if(addQTL) difCorCntProfile[,3] <- difCorCntProfile[,3]*(dtmax/qtlmax)
  invisible(difCorCntProfile)
}

plotExpressionAtMarker <- function(cross, marker, pheno.col,...){
  plotdata <- pull.pheno(cross)[,pheno.col]
  markerdata <- pull.geno(cross)[,marker]
  boxplot(plotdata[which(markerdata==1)],plotdata[which(markerdata==2)],...)
  invisible(list(plotdata[which(markerdata==1)],plotdata[which(markerdata==2)]))
}

#Plot all the profiles of traits that show differential correlation at a certain marker
plotDifCorAtMarker <- function(cross, difCntMatrix, significant=212, lodthreshold=4, marker="YBR008C_211"){
  signdifcor <- names(which(difCntMatrix[marker,] > significant))
  cat("- Starting QTL scan of",length(signdifcor),"\n")
  cross <- calc.genoprob(cross)
  res <- scanone(cross,pheno.col=signdifcor)
  aboveqtl <- names(which(apply(res[,3:ncol(res)],2,max) > lodthreshold))
  cat("- Start of",length(aboveqtl),"plots\n")
  if(length(aboveqtl) > 0){
    for(x in aboveqtl){
      plotDifCntProfile(cross,difCntMatrix,x,addQTL=TRUE)
    }
  }
  invisible(res)
}

#plot comparison of genomewide significant eQTL and eQCL
plotComparison <- function(cross, difCntMatrix, scanoneMatrix, significant=212, lodthreshold=4){
  summaryQCL <- lodscorestoscanone(cross,apply(difCntMatrix,1,function(x){sum(x > s1)}))
  summaryQTL <- lodscorestoscanone(cross,apply(scanoneMatrix[,3:ncol(scanoneMatrix)],1,function(x){sum(x > s2)}))
  plot(summaryQTL,summaryQCL, col=c("black","red"), ylab="# above threshold", main="Genomewide summary of eQTL/eQCL")
}

#Plot the expression specified by pheno.col but splits it based on the genotype at marker
plotExpressionAtMarker <- function(cross, pheno.col=1, marker="YBR008C_211"){
  genotype <- pull.geno(cross)[,marker]
  exp1 <- pull.pheno(cross)[genotype==1,pheno.col]
  exp2 <- pull.pheno(cross)[genotype==2,pheno.col]
  boxplot(exp1,exp2,main=paste("Expression of:",pheno.col),sub=paste("Split by genotype at marker:",marker))
  invisible(list(exp1,exp2))
}

#Plot the detailed correlations and difCor score of a given phenotype in a difCor object, when threshold is specified only the the elements above the threshold are plotted
plotDifCorDetail <- function(difCor, pheno.col=1, difCorThreshold=0){
  sorted <- sort(abs(difCor[[1]][,pheno.col]),decreasing=T)
  sorted <- sorted[which(sorted > difCorThreshold)]
  ordering <- names(sorted)
  plot(c(1,length(ordering)),c(0,1),xlab="Other phenotypes",ylab="Correlation",main=paste("Correlation probe ",pheno.col," at marker",attr(difCor,"marker"),sep=" "),type='n')
  points(difCor[[2]][ordering,pheno.col]^2,pch=20,col='red',type='s')
  points(difCor[[3]][ordering,pheno.col]^2,pch=20,col='green',type='s')
  points(abs(difCor[[1]][ordering,pheno.col]),pch=20,col='blue',type='l',lwd=4)
  legend("topright",c("Genotype 1","Genotype 2","Absolute difference"),pch=20,col=c("red","green","blue"))
  invisible(ordering)
}
