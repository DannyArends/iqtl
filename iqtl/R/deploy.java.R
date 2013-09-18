#
#
# deploy.java.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions to create JAVA datafile and JAVA QTLviewer
#
#

deploy.java <- function(multiresult,cross,location){
  #Copy of the 'static' scripts/svg/html + Generation of the datafile
  if(missing(multiresult)) stop("Please supply QTL scanning results from scanall")
  if(missing(cross)) stop("No Crossobject")
  iqtldatadir <- paste(installed.packages()[which(rownames(installed.packages())=="iqtl"),"LibPath"],"/iqtl/data",sep="")
  if(missing(location)){
    location <- tempdir()
    cat("- Location unavailable changing to tempdir\n")
  }else{
    cat("- Location available\n")
  }
  setwd(location)
  if(!file.exists("viewer")){
    cat("- Creating directory viewer ...\n")
    dir.create("viewer")
  }
  if(strsplit(location,"")[[1]][nchar(location)]=="/"){
    location <- paste(location,"viewer",sep="")
  }else{
    location <- paste(location,"/viewer",sep="")
  }
  setwd(iqtldatadir)
  if(file.copy("QTLviewer.jar",paste(location,"/QTLviewer.jar",sep=""))){
    cat("- Copied viewer ...\n")
  }else{
    cat("- Using old viewer ...\n")
  }
  setwd(location)
  if(file.exists("data.dat")){
    cat("- Deleting old data.dat ...\n")
    file.remove("data.dat")
  }
  datafile.java(multiresult,cross)
  cat("- data.dat generated from multiresult.\n")
  testJAVA <- system(paste("java",sep="",collapse=""),wait = TRUE,intern=T)
  if(substr(testJAVA[[1]],1,5)=="Usage"){
    cat("# Java found.\n",sep="")
    system(paste("java -jar QTLviewer.jar",sep="",collapse=""),wait = FALSE)
    cat("# Starting viewer\n",sep="")
  }else{
    cat("# Java NOT found. please install the Java Runtime Environment, and add it to your path\n",sep="")
    cat("# Please open ",location,"/QTLviewer.jar\n",sep="")
    cat("# Make sure JAVA JRE is available.\n",sep="")
    cat("#\n",sep="")
  }
  setwd(R.home())
}


datafile.java <- function(multiresult=NULL,cross=NULL){
  if(is.null(multiresult)) stop("No MQM results")
  if(is.null(cross)) stop("No Crossobject")
  mfile <- file("data.dat", "w")
  cat("individuals=",nind(cross),"\n", file = mfile,sep="")
  cat("markers=",nrow(multiresult[[1]]),"\n", file = mfile,sep="")
  cat("chromosomes=",length(unique(multiresult[[1]][,1])),"\n", file = mfile,sep="")
  l <- NULL
  slv <- 0
  for(x in unique(multiresult[[1]][,1])){
    l <- c(l,max(multiresult[[1]][which(multiresult[[1]][,1]==x),2]))
    slv <- c(slv,sum(l))
  }
  cat("length=",paste(slv,collapse=","),"\n", file = mfile,sep="")
  cat("traits=",length(multiresult),"\n", file = mfile,sep="")
  cat("minQTL=",ceiling(min(unlist(lapply(multiresult,FUN=getThird)))),"\n\n", file = mfile,sep="")  
  cat("maxQTL=",ceiling(max(unlist(lapply(multiresult,FUN=getThird)))),"\n\n", file = mfile,sep="")
  for(x in 1:length(multiresult)){
    cat("qtldata.",x-1,"=",strsplit(colnames(multiresult[[x]])[3]," ")[[1]][2],",", paste(round(multiresult[[x]][,3],digits=2),collapse=","),"\n", file = mfile,sep="")
  }
  cat("\n", file = mfile)     
  for(x in 1:nrow(multiresult[[1]])){
    cat("mapdata.",(x-1),"=",rownames(multiresult[[1]])[x],",",multiresult[[1]][x,1],",",multiresult[[1]][x,2],"\n", file = mfile,sep="")
  }
  if(!is.null(cross) && !is.null(cross$locations)){
  cat("\nlocdata=1\n", file = mfile,sep="")
  for(x in 1:length(cross$locations)){
    cat("loc.",(x-1),"=",rownames(cross$locations[[x]])[1],",",cross$locations[[x]][1,1],",",cross$locations[[x]][1,2],"\n", file = mfile,sep="")
  }
  }else{
    cat("\nlocdata=0\n", file = mfile,sep="")
  }
  if(!is.null(cross)){
    cat("\nmodeldata=1\n", file = mfile,sep="")  
    for(x in 1:length(multiresult)){
      if(!is.null(attr(multiresult[[x]],"mqmmodel"))){
        cof <- as.numeric(rownames(multiresult[[1]])%in%attr(multiresult[[x]],"mqmmodel")[[2]])
        cat("model.",(x-1),"=",paste(cof,collapse=","),"\n", file = mfile,sep="")
      }else{
        cat("model.",(x-1),"=",paste(rep(0,nrow(multiresult[[1]])),collapse=","),"\n", file = mfile,sep="")
      }
    }
  }else{
    cat("\nmodeldata=0\n", file = mfile,sep="")
  }
  close(mfile)
}

