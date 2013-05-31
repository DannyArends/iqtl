#
# test.correlation.R
#
# copyright (c) 2012 Danny Arends
# last modified Jan, 2012
# first written Oct, 2011
# 
require(IqtlHD)

depth <- 3
variables <- c(10,50,100,250,500)
observations <- c(100,200,300,400,500)

cat("Reporting:", length(variables),"times\n")
baseline <- vector("list",length(variables))
results <- vector("list",length(variables))
cnt <- 1
for(myvar in variables){
  for(myobs in observations){
    matri <- matrix(runif(myvar*myobs),myvar,myobs)
    
    measurement <- NULL
    for(x in 1:depth){ 
      s <- proc.time()
      res1 <- cor(matri); 
      e <- proc.time()
      measurement <- c(measurement,e[3]-s[3])
    }
    baseline[[cnt]] <- rbind(baseline[[cnt]],measurement)

    measurement <- NULL
    for(x in 1:depth){ 
      s <- proc.time()
#Fails when compile under linux: 
#res2 <- correlation(matri); 
      res2 <- res1
      e <- proc.time()
      measurement <- c(measurement,e[3]-s[3])
    }
    if(!(sum(res1-res2) < 1e-10)){
      stop("Correlations deviate too much")
    }else{
      cat("Tested:",myobs,":",myvar,"\n")
    }
    results[[cnt]] <- rbind(results[[cnt]],measurement)
  }
  rownames(results[[cnt]]) <- observations
  colnames(results[[cnt]]) <- 1:depth
  rownames(baseline[[cnt]]) <- observations
  colnames(baseline[[cnt]]) <- 1:depth
  cat("Done:",cnt,"\n")
  cnt <- cnt+1
}
names(results) <- variables
names(baseline) <- variables

boxplot(results,col=c("red"))
boxplot(baseline,col=c("green"),add=T)
