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
  starts <- NULL
  st  <- proc.time()
  if(saveres){
    cat("",file="snpout.txt")
    fp <- file("snpout.txt","w")
  }
  for(x in toAnalyse){
    startloc <- (x-halfregion)
    if(startloc < 1){
      left <- abs(startloc)
      startloc <- 1
    }else{
      startloc <- startloc+1
      left <- 0
    }
    endloc <- (x+halfregion+left)
    if(endloc > nrow(genotypes)){
      startloc <- startloc - (endloc - nrow(genotypes))
      endloc <- nrow(genotypes)
    }
    lodscores <- test.single.snp(as.numeric(genotypes[x,]),genotypes[startloc:endloc,])
    if(saveres) cat(x,"\t",paste(round(lodscores,1),collapse="\t"),"\n",sep="",file=fp)
    starts <- c(starts,startloc)
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

getAlleleFreq <- function(genotypes, freq = 30){
  allele_freq <- apply(geno,1,function(x){sum(as.numeric(x))}) / ncol(geno)
  invisible(allele_freq)
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

plot.snp <- function(scores, selected, map_data, chr){
  colordata <- map_data[selected,1]
  if(!missing(chr)){
    idx <- which(map_data[selected,1]==chr)
    colordata <- colordata[idx]
  }
  plot(get.chr(scores, selected, map_data, chr)[,2:3],pch=20,t='h',lwd=0.01,col=colordata)
}

get.chr <- function(scores, selected, map_data, chr){
  plotdata <- cbind(chr=map_data[selected,1], loc=map_data[selected,3], score=scores[selected])
  rownames(plotdata) <- rownames(map_data[selected,])
  if(!missing(chr)){
    idx <- which(map_data[selected,1]==chr)
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
    cat(x,nrow(genotypes),region, locations,"\n")
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
  mm
}

plot_test <- function(genotypes, pcors, myrange){
  if(missing(myrange)){
    subselection <- 1:nrow(genotypes)
  }else{
    subselection <- which(genotypes[,2] < myrange[2] & genotypes[,2] > myrange[1])
  }
  pcors <- pcors[subselection,]
  bprange <- c(min(genotypes[subselection,2]),max(genotypes[subselection,2]))
  plot(bprange,bprange,t='n')
  region <- attr(pcors,"region")
  halfregion <- region/2
  heat11 <- rgb(seq(0,1,0.1),seq(1,0,-0.1),seq(1,0,-0.1),seq(0,1,0.1))
  colorz <- apply((round(10*abs(pcors),0)+1),2,function(x){heat11[as.numeric(x)]})
  cat(dim(colorz),"\n")
  for(x in subselection){
    locations <- getLocs(x, halfregion, nrow(genotypes))
    points(rep(genotypes[x,2],region),genotypes[locations[1]:locations[2],2],pch=20,col=colorz[x,],cex=0.1)
  }
}

setwd("E:/GBIC/LFN")
load("genotypes.Rdata")
load("snp.test.Rdata")
map_data <- read.table("map_lfn_005.txt",sep=",",row.names=1,header=TRUE)
#result <- snp.test(genotypes, region=50, update.time=100)

scores   <- get.above(result, 100, 5)
selected <- which(!is.na(scores))

#scores <- get.above(result, 3, 0)
#selected <- which(map_data[,1]==1)

plot.snp(scores, selected, map_data)

my_chr1 <- get.genotypes(genotypes, selected, map_data, 1)
pcors   <- get.pcorrelation(my_chr1, region=1500)

jpeg("chr1_new2mb.jpg")
plot_test(my_chr1,pcors,c(1,2000000))
dev.off()

