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
  cat("",file="snpout.txt")
  if(saveres) fp <- file("snpout.txt","w")
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
    if(endloc > nrow(genotypes)) endloc <- nrow(genotypes)
    lodscores <- test.single.snp(as.numeric(genotypes[x,]),genotypes[startloc:endloc,])
    cat(x,"\t",paste(round(lodscores,1),collapse="\t"),"\n",sep="",file=fp)
    starts <- c(starts,startloc)
    if(verbose && x %% update.time == 0){
      cat("Analysis of",paste("(",x-update.time,"..",x,")",sep=""),"in:",(proc.time()-st)[3],"seconds\n")
      st  <- proc.time()
    }
  }
  close(fp)
  mm <- read.table("snpout.txt",sep="\t",row.names=1)
  attr(mm,"region") = region
  mm
}

getAlleleFreq <- function(genotypes, freq = 30){
  allele_freq <- apply(geno,1,function(x){sum(as.numeric(x))}) / ncol(geno)
  invisible(allele_freq)
}

get.above <- function(result, lod.cutoff = 100){
  apply(apply(result,2,as.numeric) > lod.cutoff,1,sum)
}

removeInf <- function(res){
  apply(res,2,function(x){x[which(x==Inf)] <- NA;x})
}

get.above.raw <- function(result){
  result <- removeInf(result)
  apply(apply(result,2,as.numeric),1,sum,na.rm=TRUE)
}

#map_data <- read.table("map_lfn_005.txt",sep=",",row.names=1,header=TRUE)

res <- snp.test(genotypes, region=50, update.time=100)

get.above(res)

plot(res[[1]]/res[[2]],t='p',col=((res[[1]]/res[[2]]) < 0.8)+1, pch=20)

goodsnp <- which(((res[[1]]/res[[2]]) > 0.8))

goodones <- which(get.above.raw(removeInf(res)) > 400)


locs <- map_data[rownames(genotypes[goodones,]),1:2]

plot(cbind(x=locs[,2],y=get.above.raw(removeInf(res[goodones,]))),pch=20,t='h')

plot(cbind(x=locs[,2],y=rep(1,nrow(locs))),pch=20)

oldgeno <- geno[1:1000,]
newgeno <- geno[goodsnp,]
