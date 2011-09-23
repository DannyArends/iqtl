#
#
# deplot.svg.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions to create SVG datafile and SVG QTLviewer
#
#


deploy.svg <- function(multiresult,cross,location,FireFox="C:/Program Files/Mozilla Firefox/firefox.exe"){
  #Copy of the 'static' scripts/svg/html + Generation of the datafile
  if(missing(multiresult)) stop("Please supply QTL scanning results from scanall")
  iqtldatadir <- paste(installed.packages()[which(rownames(installed.packages())=="iqtl"),"LibPath"],"/iqtl/data",sep="")
  if(missing(location)){
    location <- tempdir()
    cat("- Location unavailable changing to homedir\n")
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
  if(file.copy("functions.js",paste(location,"/functions.js",sep=""))){
    cat("- Copied functions.js ...\n")
  }else{
    cat("- functions.js already exists. Not overwriting.\n")
  }
  if(file.copy("qtlviewer.svg",paste(location,"/qtlviewer.svg",sep=""))){
    cat("- Copied qtlviewer.svg ...\n")
  }else{
    cat("- qtlviewer.svg already exists. Not overwriting.\n")
  }
  if(file.copy("index.html",paste(location,"/index.html",sep=""))){
    cat("- Copied index.html ...\n")
  }else{
    cat("- index.html already exists. Not overwriting.\n")
  }
  setwd(location)
  if(file.exists("data.js")){
    cat("- Deleting old data.js ...\n")
    file.remove("data.js")
  }
  datafile.svg(multiresult,cross)
  cat("- data.js generated from multiresult.\n")
  if(!is.null(FireFox) && file.exists(FireFox)){
    system(paste("\"",FireFox,"\" file://",location,"/index.html",sep="",collapse=""),wait = FALSE)
    cat("# Starting FireFox\n",sep="")
  }else{
    cat("#\n",sep="")
    cat("# Please open ",location,"/index.html using Mozilla FireFox.\n",sep="")
    cat("# Any webbrowser capable of handling SVG 1.1 can be used.\n",sep="")
    cat("#\n",sep="")
  }
  setwd(R.home())
}

datafile.svg <- function(multiresult=NULL,cross=NULL){
  #Generation of the datafile
  if(is.null(multiresult)) stop("Please supply QTL scanning results from scanall")
  mfile <- file("data.js", "w")
  cat("//\n//\n// data.js\n//\n// copyright (c) 2010, Danny Arends\n// last modified mrt, 2010\n// first written mrt, 2010\n//\n// Generated Datafile for QTLviewer\n// Please do NOT edit...\n//\n//\n", file = mfile,sep="")
  cat("var qtldata=[];\nvar mapdata=[];\nvar chrL=[];\n", file = mfile,sep="")
  cat("var markers=",nrow(multiresult[[1]]),";\n", file = mfile,sep="")
  cat("var chromosomes=",length(unique(multiresult[[1]][,1])),";\n", file = mfile,sep="")
  l <- NULL
  slv <- 0
  for(x in unique(multiresult[[1]][,1])){
    l <- c(l,max(multiresult[[1]][which(multiresult[[1]][,1]==x),2]))
    slv <- c(slv,sum(l))
  }
  cat("var chrL=[",paste(slv,collapse=","),"];\n", file = mfile,sep="")
  cat("var traits=",length(multiresult),";\n", file = mfile,sep="")
  cat("var maxQTL=",max(unlist(lapply(multiresult,FUN=getThird))),";\n\n", file = mfile,sep="")
  for(x in 1:length(multiresult)){
    cat("qtldata[",(x-1),"] = [\"",strsplit(colnames(multiresult[[x]])[3]," ")[[1]][2],"\",", paste(round(multiresult[[x]][,3],dig=2),collapse=","),"];\n", file = mfile,sep="")
  }
  cat("\n", file = mfile)     
  for(x in 1:nrow(multiresult[[1]])){
    cat("mapdata[",(x-1),"] = [\"",rownames(multiresult[[1]])[x],"\",",multiresult[[1]][x,1],",",multiresult[[1]][x,2],"];\n", file = mfile,sep="")
  }
  if(!is.null(cross) && !is.null(cross$locations)){
    cat("\nvar locdata=[];\n", file = mfile,sep="")
    for(x in 1:length(cross$locations)){
      cat("locdata[",(x-1),"] = [\"",rownames(cross$locations[[x]])[1],"\",",cross$locations[[x]][1,1],",",cross$locations[[x]][1,2],"];\n", file = mfile,sep="")
    }
  }else{
    cat("\nvar locdata=null;\n", file = mfile,sep="")
  }
  if(!is.null(cross)){
    cat("\nvar modeldata=[];\n", file = mfile,sep="")  
    for(x in 1:length(multiresult)){
      if(!is.null(attr(multiresult[[x]],"mqmmodel"))){
        cof <- as.numeric(rownames(multiresult[[1]])%in%attr(multiresult[[x]],"mqmmodel")[[2]])
        cat("modeldata[",(x-1),"] = [",paste(cof,collapse=","),"];\n", file = mfile,sep="")
      }else{
        cat("modeldata[",(x-1),"] = [",paste(rep(0,nrow(multiresult[[1]])),collapse=","),"];\n", file = mfile,sep="")
      }
    }
  }else{
    cat("\nvar modeldata=null;\n", file = mfile,sep="")
  }
  close(mfile)
}

