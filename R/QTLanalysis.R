#####################################################################
#
# QTLanalysis.R
#
# copyright (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Feb, 2009
# 
# Part of the ClusterJobs package
# Contains: run_cluster
#
######################################################################

######################################################################
#
# generateQTLfile: Generates a QTL file which is executed (R CMD BATCH <qtlfile>) by a cluster
# generateESTfile: Generates a QTL file to estimate runtime for each item
#
######################################################################

generateQTLfile <- function(dbpath = "", job, b_size,Ntraits,name,jobid,map,method,model,step,investigationname=""){
	#Generate a qtl<#>.R file to do the analysis on
	qtlfile <- paste("~/run",jobid,"/subjob",job,".R",sep="")
	#load needed libraries
	cat("library(qtl,lib.loc='~/libs')","\n",sep="",file=qtlfile)
	cat("library(bitops,lib.loc='~/libs')","\n",sep="",file=qtlfile,append=T)
	cat("library(RCurl,lib.loc='~/libs')","\n",sep="",file=qtlfile,append=T)
    cat("library(ClusterJobs,lib.loc='~/libs')","\n",sep="",file=qtlfile,append=T)
    cat("source(\"",paste(dbpath,"/api/R",sep=""),"\")\n",sep="",file=qtlfile,append=T)
	#Print our report function
	cat("\nreport <- function(status,text){\n",file=qtlfile,append=T)
	cat("\ttask <- ",jobid,"\n",file=qtlfile,append=T)
	cat("\ttext <- substr(URLencode(text),0,100)\n",file=qtlfile,append=T)
	cat("\tjob <- ",job,"\n",file=qtlfile,append=T)
	cat("\tlink <- paste(\"",dbpath,"/taskreporter?job=\",task,\"&subjob=\",job,\"&statuscode=\",status,\"&statustext=\",text,sep=\"\")\n",sep="",file=qtlfile,append=T)
	cat("\tgetURL(link)\n",file=qtlfile,append=T)
	cat("\tif(status==-1){\n\t\tcat(\"!!!\",text,\"!!!\")\n\t\tq(\"no\")\n\t}\n",file=qtlfile,append=T)
	cat("}\n\n",file=qtlfile,append=T)
	#Start by sending a message (so we know we're running
	cat("report(2,\"LoadingCrossobject\")\n",file=qtlfile,append=T)
	cat(Generate_Statement(paste("load(\"~/run",jobid,"/cross.RData\")","\n",sep="")),file=qtlfile,append=T)
	cat("report(2,\"FinishedLoading\")\n",file=qtlfile,append=T)
    if(step!=0){
    cat(Generate_Statement(paste("cross <- calc.genoprob(cross, step=",step,")","\n",sep="")),file=qtlfile,append=T)
  }
  cat("report(2,\"FinishedMarkerImputation\")\n",file=qtlfile,append=T)
	#Which genotypes do we actually need ? (throw away any remaining)
	if(((job-1)*b_size)+b_size > Ntraits){
		#final batch
		cat(Generate_Statement(paste("cross$pheno <- cross$pheno[",((job-1)*b_size)+1,":",Ntraits,"]","\n",sep="")),file=qtlfile,append=T)
	}else{
		cat(Generate_Statement(paste("cross$pheno <- cross$pheno[",((job-1)*b_size)+1,":",((job-1)*b_size)+b_size,"]","\n",sep="")),file=qtlfile,append=T)
	}
	#Calculate results
  #this is a hack
  #TODO request map method paramters and 
  if(map != "cimall"){
    cat(Generate_Statement(paste("results <- ",map,"(cross,method='",method,"',model='",model,"',multicore=FALSE,verbose=TRUE)","\n",sep="")),file=qtlfile,append=T)
  }else{
    cat(Generate_Statement(paste("results <- ",map,"(cross,method='",method,"',multicore=FALSE,verbose=TRUE)","\n",sep="")),file=qtlfile,append=T)
  }
	cat("report(2,\"StoringQTLresults\")\n",file=qtlfile,append=T)
  #Store Results
	cat(Generate_Statement(paste("ResultsToMolgenis(result=results,resultname='",name,"',Trait_num=",(job-1)*b_size,",investigationname=\"",investigationname,"\")","\n",sep="")),file=qtlfile,append=T)
	cat("report(3,\"JobFinished\")\n",file=qtlfile,append=T)
	#Quit
	cat("q(\"no\")","\n",sep="",file=qtlfile,append=T)
}

generateESTfile <- function(dbpath = "", map, method, model, step, jobid){
	#Generates a qtlfile to estimate walltime
	qtlfile <- paste("~/run",jobid,"/ESTtime.R",sep="")
	#Print our report function
	cat("\nreport <- function(status,text){\n",file=qtlfile)
	cat("\ttask <- ",jobid,"\n",file=qtlfile,append=T)
	cat("\ttext <- substr(URLencode(text),0,100)\n",file=qtlfile,append=T)
	cat("\tlink <- paste(\"",dbpath,"/taskreporter?job=\",task,\"&subjob=0&statuscode=\",status,\"&statustext=\",text,sep=\"\")\n",sep="",file=qtlfile,append=T)
	cat("\tgetURL(link)\n",file=qtlfile,append=T)
	cat("\tif(status==-1){\n\t\tcat(\"!!!\",text,\"!!!\")\n\t\t\n\t\tq(\"no\")\n\t}\n",file=qtlfile,append=T)
	cat("}\n\n",file=qtlfile,append=T)
	#load needed libraries
	cat("library(qtl,lib.loc='~/libs')","\n",sep="",file=qtlfile,append=T)
	cat("library(bitops,lib.loc='~/libs')","\n",sep="",file=qtlfile,append=T)	
	cat("library(RCurl,lib.loc='~/libs')","\n",sep="",file=qtlfile,append=T)
	cat(Generate_Statement(paste("load(\"~/run",jobid,"/cross.RData\")","\n",sep="")),file=qtlfile,append=T)
	cat("cross$pheno <- cross$pheno[1]","\n",sep="",file=qtlfile,append=T)
  if(step!=0){
    cat(Generate_Statement(paste("cross <- calc.genoprob(cross, step=",step,")","\n",sep="")),file=qtlfile,append=T)
  }
  cat("report(2,\"FinishedMarkerImputation\")\n",file=qtlfile,append=T)
  if(map != "cimall"){
	  cat(Generate_Statement(paste("results <- ",map,"(cross,method='",method,"',model='",model,"',multicore=FALSE,verbose=TRUE)","\n",sep="")),file=qtlfile,append=T)
  }else{
    cat(Generate_Statement(paste("results <- ",map,"(cross,method='",method,"',multicore=FALSE,verbose=TRUE)","\n",sep="")),file=qtlfile,append=T)
  }
	#Quit
	cat("q(\"no\")","\n",sep="",file=qtlfile,append=T)
}




