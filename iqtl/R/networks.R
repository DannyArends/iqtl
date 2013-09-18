######################################################################
#
# networks.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2011
# first written mrt, 2011
# 
# R functions: mqmmodelsasnetwork, scanonetosif
#
#
######################################################################

mqmmodelsasnetwork <- function(cross,result){
  if(is.null(cross)){
		stop("No cross object. Please supply a valid cross object.") 
	}
  if(!any(class(result)=="mqmmulti")){
    stop("No mqmmulti object. Please supply a valid mqmmulti object.") 
  }
  models <- lapply(FUN=mqmgetmodel,result)
  namez <- colnames(pull.pheno(cross))
  cat(file="QTLnetwork.sif","",append=F)
  cat(file="QTLnodes.sif","",append=F)
  for(x in 1:length(models)){
    cat(file="QTLnodes.sif",namez[x],"Trait\t0\n",sep="\t",append=TRUE)
    for(y in 1:length(models[[x]][[2]])){
      cat(file="QTLnodes.sif",models[[x]][[2]][y],"Marker",models[[x]][[4]][y],"\n",sep="\t",append=TRUE)
      cat(file="QTLnetwork.sif",namez[x],"QTLeffect",models[[x]][[2]][y],result[[x]][models[[x]][[2]][y],3],"\n",sep="\t",append=TRUE)
    }
  }
}


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