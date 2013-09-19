#
# snp.test.R
#
# Copyright (c) 2012 Danny Arends
# Last modified Mar, 2012
# First written Mar, 2012
# 

test.single.snp <- function(genotype, genotypes){
  apply(genotypes,1,function(x){
    -log10(anova(lm(genotype ~ as.factor(x)))[[5]][1])
  })
}

snp.test <- function(genotypes, toAnalyse=1:nrow(genotypes), region = 50, update.time=10, saveres = TRUE, verbose = TRUE){
  halfregion <- region/2
  st  <- proc.time()
  if(saveres){
    cat("",file="snpout.txt")
    fp <- file("snpout.txt","w")
  }
  for(x in toAnalyse){
    location <- getLocs(x, halfregion, nrow(genotypes))
    lodscores <- test.single.snp(as.numeric(genotypes[x,]),genotypes[location[1]:location[2],])
    if(saveres) cat(x,"\t",paste(round(lodscores,1),collapse="\t"),"\n",sep="",file=fp)
    if(verbose && x %% update.time == 0){
      cat("Analysis of",paste("(",x-update.time,"..",x,")",sep=""),"in:",(proc.time()-st)[3],"seconds\n")
      st  <- proc.time()
    }
  }
  if(saveres) close(fp)
  mm <- read.table("snpout.txt",sep="\t",row.names=1)
  attr(mm,"region") = region
  mm
}

# Note: Markers are columns of the genotype matrix
getAlleleFreq <- function(genotypes, MARGIN = 2){
  invisible(apply(genotypes, MAGIN=MARGIN, function(x){
  	nx <- as.numeric(x)
  	sum(nx, na.rm=TRUE) / length(which(!is.na(nx)))
  }))
}

get.above <- function(result, lod.cutoff = 100, minimum = 5, zero.na = TRUE){
  cleaned <- (apply(apply(result,2,as.numeric) > lod.cutoff,1,sum)-1)
  if(zero.na) cleaned[which(cleaned < minimum)] <- NA
  invisible(cleaned)
}

removeInf <- function(res){
  apply(res,2,function(x){x[which(x==Inf)] <- NA;x})
}

get.above.raw <- function(result){
  result <- removeInf(result)
  apply(apply(result,2,as.numeric),1,sum,na.rm=TRUE)
}

plot.overview <- function(x, ...){
  args <- list(...)
  colordata <- args$map[args$selected, 1]
  if(!missing(args$chr)){
    idx <- which(args$map[args$selected, 1] == args$chr)
    colordata <- colordata[idx]
  }
  plot(get.chr(x, args$selected, args$map, args$chr)[, 2:3], pch=20, t='h', lwd=0.01, col=colordata)
}

get.chr <- function(scores, selected, map, chr){
  plotdata <- cbind(chr=map[selected, 1], loc=map[selected, 3], score=scores[selected])
  rownames(plotdata) <- rownames(map[selected,])
  if(!missing(chr)){
    idx <- which(map[selected, 1] == chr)
    plotdata <- plotdata[idx,]
  }
  invisible(plotdata)
}

get.genotypes <- function(genotypes, selected, map_data, chr){
  plotdata <- cbind(chr=map_data[selected,1], loc=map_data[selected,3], genotypes[selected,])
  plotdata <- apply(plotdata,2,as.numeric)
  rownames(plotdata) <- rownames(map_data[selected,])
  colnames(plotdata) <- c("chr","loc", colnames(genotypes))
  if(!missing(chr)){
    idx <- which(map_data[selected,1]==chr)
    plotdata <- plotdata[idx,]
    attr(plotdata,"chr") <- chr
  }
  invisible(plotdata)
}

getLocs <- function(x, halfregion, maximum){
  startloc <- (x-halfregion)
  if(startloc < 1){
    left <- abs(startloc)
    startloc <- 1
  }else{
    startloc <- startloc+1
    left <- 0
  }
  endloc <- (x+halfregion+left)
  if(endloc > maximum){
    startloc <- startloc - (endloc - maximum)
    endloc <- maximum
  }
  return(c(startloc,endloc))
}

get.pcorrelation <- function(genotypes, toAnalyse=1:nrow(genotypes), region=100, update.time=1000, saveres = TRUE, verbose = TRUE){
  st  <- proc.time()
  halfregion <- region/2
  mm <- NULL
  if(saveres){
    cat("",file="snpcorout.txt")
    fp <- file("snpcorout.txt","w")
  }
  for(x in toAnalyse){
    locations <- getLocs(x, round(region/2,0),nrow(genotypes))
    #cat(x,nrow(genotypes),region, locations,"\n")
    if(verbose && x %% update.time == 0){
      cat("Analysis of",paste("(",x-update.time,"..",x,")",sep=""),"in:",(proc.time()-st)[3],"seconds\n")
      st  <- proc.time()
    }
    corscores <- cor(genotypes[x,3:ncol(genotypes)], t(genotypes[locations[1]:locations[2],3:ncol(genotypes)]))
    if(saveres) cat(rownames(genotypes)[x],"\t",paste(round(corscores,4),collapse="\t"),"\n",sep="",file=fp)
  }
  if(saveres) close(fp)
  mm <- read.table("snpcorout.txt",sep="\t",row.names=1)
  attr(mm,"region") = region
  attr(mm,"chr") = attr(genotypes,"chr")
  mm
}

plot.snpplot <- function(x, pcors, myrange, ...){
  if(missing(myrange)){
    subselection <- 1:nrow(x)
  }else{
    subselection <- which(x[,2] < myrange[2] & x[,2] > myrange[1])
  }
  pcors <- pcors[subselection,]
  bprange <- c(min(x[subselection,2]),max(x[subselection,2]))
  cat(bprange,"\n")
  region <- attr(pcors,"region")
  chr <- attr(pcors,"chr")
  
  plot(bprange,bprange,main=paste("SNP Correlation Chr ",chr,sep=""),sub=paste("Using a region of ",region," snps",sep=""),t='n')
  halfregion <- region/2
  heat11 <- rgb(seq(0,1,0.1),seq(1,0,-0.1),seq(1,0,-0.1),seq(0,1,0.1))
  colorz <- apply((round(10*abs(pcors),0)+1),2,function(x){heat11[as.numeric(x)]})
  cat(dim(colorz),"\n")
  for(x in subselection){
    locations <- getLocs(x, halfregion, nrow(x))
    points(rep(x[x,2],region),x[locations[1]:locations[2],2],pch=20,col=colorz[x,],cex=0.1)
  }
}

plot.snps <- function(x, map, selected, chr, region, ...){
  args <- list(...)  
  genoFMT <- get.genotypes(x, selected, map, chr)

  cat("Selected",nrow(my_chr),"for plotting\n")
  pcors   <- get.pcorrelation(genoFMT, region=region)

  png(paste("Chr", chr, "_top25k.jpg", sep=""),width = 1024, height = 1024)
  plot.snpplot(genoFMT, pcors)
  dev.off()
}

#test.snp <- function(){
#  memory.limit(3000)
#  setwd("E:/GBIC/LFN")
#  load("genotypes.Rdata")
#  load("snp.test.Rdata")
#  map_data <- read.table("map_lfn_005.txt",sep=",",row.names=1,header=TRUE)
#  #result <- snp.test(genotypes, region=50, update.time=100)

#  scores   <- get.above.raw(result)
#  names(scores) <- paste("m",1:nrow(genotypes),sep="")
#  top25k <- names(sort(scores,decreasing=TRUE)[1:25000])
#  selected <- which(names(scores) %in% top25k)
#  cat("Selected",length(selected),"from",length(scores),"\n")

#  png(paste("snp_overview.jpg",sep=""),width = 800, height = 600)
#  plot.overview(scores, selected, map_data)
#  dev.off()

  #plot.snps(genotypes, selected, map_data, 1,1000)
  #plot.snps(genotypes, selected, map_data, 2,1000)
  #plot.snps(genotypes, selected, map_data, 3,1000)
  #plot.snps(genotypes, selected, map_data, 4,1000)
  #plot.snps(genotypes, selected, map_data, 5,1000)

#  my_chr <- get.genotypes(genotypes, selected, map_data)
#  FT <- read.table("FT.txt")

#  res1 <- test.single.snp(FT[,1],my_chr[,3:ncol(my_chr)])
#  res2 <- test.single.snp(FT[,2],my_chr[,3:ncol(my_chr)])
#  res3 <- test.single.snp(FT[,3],my_chr[,3:ncol(my_chr)])
#  res4 <- test.single.snp(FT[,4],my_chr[,3:ncol(my_chr)])
#  FTmapping <- rbind(res1,res2,res3,res4)
#  rownames(FTmapping) <- c("FT_r1","FT_r2","FT_r3","FT_mean")

#  png("FT_r1.jpg",width = 1024, height = 600)
#    plot(x=map_data[selected,3],y=res1,col=map_data[selected,1],pch=20,t="o",main="QTLprofile of FT_r1")
#  dev.off()

#  png("FT_r2.jpg",width = 1024, height = 600)
#    plot(x=map_data[selected,3],y=res2,col=map_data[selected,1],pch=20,t="o",main="QTLprofile of FT_r2")
#  dev.off()

#  png("FT_r3.jpg",width = 1024, height = 600)
#    plot(x=map_data[selected,3],y=res3,col=map_data[selected,1],pch=20,t="o",main="QTLprofile of FT_r3")
#  dev.off()

#  png("FT_mean.jpg",width = 1024, height = 600)
#    plot(x=map_data[selected,3],y=res4,col=map_data[selected,1],pch=20,t="o",main="QTLprofile of FT_mean")
#  dev.off()

#  RJ <- read.table("RJ.txt",sep="\t")
#  RJ <- RJ[rownames(FT),]

#  cat("",file="seed_qtls.txt")
#  fp <- file("seed_qtls.txt","w")
#  for(x in seq(1,ncol(RJ),2)){
#    phenotype  <- apply(RJ[,x:(x+1)],1,mean,na.rm=T)
#    qtl_result <- test.single.snp(phenotype, my_chr[,3:ncol(my_chr)])
#    cat(substring(colnames(RJ)[x], 0, nchar(colnames(RJ)[x])-2),"\t",paste(round(qtl_result,3),collapse="\t"),"\n",sep="",file=fp)
#    #png(paste("plotRJ",x,".jpg",sep=""),width = 1024, height = 600)
#    plot(x=map_data[selected,3],y=qtl_result,col=map_data[selected,1],pch=20,t="o",main=paste("QTLprofile of mean (",paste(colnames(RJ)[x:(x+1)],collapse=", "),")"))
#    #dev.off()
#  }
#  close(fp)

#  genes = apply(is.na(map_data[,4:5]),1,sum) + 1
#  points(map_data[which(map_data[,3]>=1),2],-1000*genes,pch=20,cex=0.05,col=genes)
#  dev.off()
#}
