#
# backwardelimination.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: mqmmultitoscanone, mqmmodelsasnetwork
#  Basic scripts for Causal inference
#

inverseF <- function(df1,df2,alpha=0.05){
  m <- NULL
  for(x in 1:length(df1)){
    r <- NULL
    for(y in 1:length(df2)){
    result <- .C("inverseF_R",df1=as.integer(df1[x]),
                              df2=as.integer(df2[y]),
                              alpha=alpha,
                              out=0)
    r <- c(r,result$out)
    }
    m <- rbind(m,r)
  }
  rownames(m) <- df1
  colnames(m) <- df2
  m
}

backwardeliminate <- function(cross,cofactors,alpha=0.05){
  
}