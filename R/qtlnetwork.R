#
# qtlnetwork.R
#
# copyright (c) 2010, Danny Arends
# last modified jan, 2011
# first written jan, 2011
# 
# R functions: scanonetosif
#
# Basic scripts creating networks from scanone results
#

scanonetosif <- function(results, lodcutoff=3.0, verbose=TRUE){
  date <- paste(format(Sys.Date(), "%d%b%Y"),"-",format(Sys.time(), "%H-%M"),sep="")
  cat(file=paste("qtlnetwork",date,".nodes",sep=""),"",append=F)
  if(verbose) cat(date,"Creating: qtlnetwork.nodes\n")
  for(marker in 1:nrow(results)){
    cat(file=paste("qtlnetwork",date,".nodes",sep=""),rownames(results)[marker],"\t",results[marker,1],"\t",results[marker,2],"\tMARKER\n",append=TRUE,sep="")  
  }
  for(trait in 3:ncol(results)){
    cat(file=paste("qtlnetwork",date,".nodes",sep=""),colnames(results)[trait],"\t0\t0\tTRAIT\n",append=TRUE,sep="") 
  }
  if(verbose) cat(date,"Creating: qtlnetwork.sif\n")
  cat(file=paste("qtlnetwork_",date,".sif",sep=""),"",append=F)
  cat(file=paste("qtlnetwork",date,".edge",sep=""),"",append=F)
  for(trait in 3:ncol(results)){
    for(marker in 1:nrow(results)){
     if(results[marker,trait] > lodcutoff){
        cat(file=paste("qtlnetwork_",date,".sif",sep=""),rownames(results)[marker],"\t",paste(marker,"-",trait,sep=""),"\t",colnames(results)[trait],"\n",append=TRUE,sep="")
        cat(file=paste("qtlnetwork",date,".edge",sep=""),paste(rownames(results)[marker]," (",paste(marker,"-",trait,sep=""),") ",colnames(results)[trait],sep=""),"\t",results[marker,trait],"\n",append=TRUE,sep="")
      }
    }
  }
  if(verbose) cat(date,"Done open the qtlnetwork.sif file in Cytoscape (r)\n")
}

