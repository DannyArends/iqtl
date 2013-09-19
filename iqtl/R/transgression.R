#
# transgression.R
#
# copyright (c) 2012, Danny Arends
# last modified Apr, 2012
# first written Apr, 2012
# 
# R functions: plot.transgression, test.plot.transgression
#
plot.transgression <- function(x, labels = c("Bay","Sha"), col = c("red","blue"), reorder = TRUE, ...){
  bayidx <- grep(labels[1], colnames(x))
  shaidx <- grep(labels[2], colnames(x))
  
  #Extract the parental expression data from the mdata object
  rils <- x[,-c(bayidx,shaidx)]
  bay <- x[,bayidx]
  sha <- x[,shaidx]

  baymax <- apply(bay,1,max,na.rm=T)
  baymin <- apply(bay,1,min,na.rm=T)
  shamax <- apply(sha,1,max,na.rm=T)
  shamin <- apply(sha,1,min,na.rm=T)

  minima <- apply(cbind(baymin,shamin),1,min)
  maxima <- apply(cbind(baymax,shamax),1,max)
  idx <- 1
  sizes <- apply(rils,1,function(x){
    r <- length(which(x < minima[idx]))+length(which(x > maxima[idx]))
    idx<<-idx+1
    r
  })
  if(reorder){
    #Re-ordering, showing the traits with most transgression on top
    norder <- names(sort(sizes))

    bay  <- bay[norder,]
    sha  <- sha[norder,]
    rils <- rils[norder,]

    baymax <- apply(bay,1,max,na.rm=T)
    baymin <- apply(bay,1,min,na.rm=T)
    shamax <- apply(sha,1,max,na.rm=T)
    shamin <- apply(sha,1,min,na.rm=T)
  }
  main = "Transgression in RIL"
  ylab = "Metabolite ID"
  xlab = "Scaled expression [-1, 1]"
  plot(c(-5 , 5), c(0, nrow(rils)), t='n', main=main, ylab=ylab, xlab=xlab)
  for(x in 1:nrow(rils)){
    minn <- min(baymin[x],shamin[x])
    maxx <- max(baymax[x],shamax[x])
    mrange <- 0.5*(maxx-minn)
    points(y=rep(x,ncol(rils)),x=((rils[x,]-(minn))/mrange)-1,col=rgb(0,0,0,0.5),pch=16,cex=0.1)
    points(y=c(x,x),x=((c(baymin[x],baymax[x])-(minn))/mrange)-1,col=col[1],pch=16,cex=0.7)
    points(y=c(x,x),x=((c(shamin[x],shamax[x])-(minn))/mrange)-1,col=col[2],pch=16,cex=0.7)
  }
  legend("bottomleft",c("Bay-0","Sha","RIL"),col=c(col,"black"),pch=16)
  abline(v=c(-1,1),lty=2)
}

test.plot.transgression <- function(){
  setwd("E:\\GBIC\\Ronny Joosen")
  plot.transgression(read.table("transgression.txt"))
}
