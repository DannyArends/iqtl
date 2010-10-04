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
    result <- .C("lodscorebyem_R",df1=as.integer(df1),
                                df2=as.integer(df2),
                                alpha=alpha,
                                out=0)
  result
}

backwardeliminate <- function(cross,cofactors,alpha=0.05){

}