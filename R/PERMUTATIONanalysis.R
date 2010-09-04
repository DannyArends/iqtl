#####################################################################
#
# PERMUTATIONanalysis.R
#
# copyright (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Nov, 2009
# 
# Part of the ClusterJobs package
# Contains: run_cluster
#
######################################################################

######################################################################
#
# generatePERMfile: Generates a QTL file which is executed (R CMD BATCH <qtlfile>) by a cluster
# generateESTfile: Generates a QTL file to estimate runtime for each item
#
######################################################################

generatePERMfile <- function(DBpath = "", job, b_size,Ntraits,name,taskID,nperm,map,method,model){
	#Generate a qtl<#>.R file to do the analysis on
	qtlfile <- paste("~/run",taskID,"/subjob",job,".R",sep="")
	#load needed libraries
	cat("library(qtl,lib.loc='~/libs')","\n",sep="",file=qtlfile)
	cat("library(bitops,lib.loc='~/libs')","\n",sep="",file=qtlfile,append=T)
	cat("library(RCurl,lib.loc='~/libs')","\n",sep="",file=qtlfile,append=T)
	cat("library(snow,lib.loc='~/libs')","\n",sep="",file=qtlfile,append=T)
cat("source(\"",paste(DBpath,"/api/R",sep=""),"\")\n",file=qtlfile,append=T)
	#Print our report function
	cat("\nreport <- function(status,text){\n",file=qtlfile,append=T)
	cat("\ttask <- ",taskID,"\n",file=qtlfile,append=T)
	cat("\ttext <- substr(URLencode(text),0,100)\n",file=qtlfile,append=T)
	cat("\tjob <- ",job,"\n",file=qtlfile,append=T)
	cat("\tlink <- paste(\"",DBpath,"/taskreporter?job=\",task,\"&subjob=\",job,\"&statuscode=\",status,\"&statustext=\",text,sep=\"\")\n",sep="",file=qtlfile,append=T)
	cat("\tgetURL(link)\n",file=qtlfile,append=T)
	cat("\tif(status==-1){\n\t\tcat(\"!!!\",text,\"!!!\")\n\t\tq(\"no\")\n\t}\n",file=qtlfile,append=T)
	cat("}\n\n",file=qtlfile,append=T)
	#Start by sending a message (so we know we're running
	cat("report(2,\"LoadingCrossobject\")\n",file=qtlfile,append=T)
	cat(Generate_Statement(paste("load(\"~/run",taskID,"/cross.RData\")","\n",sep="")),file=qtlfile,append=T)
	cat("report(2,\"FinishedLoading\")\n",file=qtlfile,append=T)
	#Which genotypes do we actually need ? (throw away any remaining)
	if(((job-1)*b_size)+b_size > Ntraits){
		#final batch
		cat(Generate_Statement(paste("cross$pheno <- cross$pheno[",((job-1)*b_size)+1,":",Ntraits,"]","\n",sep="")),file=qtlfile,append=T)
	}else{
		cat(Generate_Statement(paste("cross$pheno <- cross$pheno[",((job-1)*b_size)+1,":",((job-1)*b_size)+b_size,"]","\n",sep="")),file=qtlfile,append=T)
	}
	#Store results
	cat("Thresholds <- NULL\n",file=qtlfile,append=T)
	cat("Names <- NULL\n",file=qtlfile,append=T)
	cat("for(num in 1:nphe(cross)){\n",file=qtlfile,append=T)
	cat(Generate_Statement(paste("bresults <- bootstrap(cross,n.run=",nperm,",pheno.col=num,Funktie=",map,",method='",method,"',model='",model,"',multicore=FALSE,verbose=TRUE)","\n",sep="")),file=qtlfile,append=T)
	cat(Generate_Statement(paste("results <- summary(MQMpermObject(bresults))[,1]","\n",sep="")),file=qtlfile,append=T)
	cat(Generate_Statement(paste("Thresholds <- rbind(Thresholds,results)","\n",sep="")),file=qtlfile,append=T)
	cat(Generate_Statement(paste("Names <- c(Names,substr(colnames(bresults[[1]])[3],5,nchar(colnames(bresults[[1]])[3])))","\n",sep="")),file=qtlfile,append=T)
	cat("}\n",file=qtlfile,append=T)
	cat("rownames(Thresholds) <- Names\n",file=qtlfile,append=T)
	cat("Thresholds\n",file=qtlfile,append=T)
	cat(Generate_Statement(paste("PermToMolgenis(permResults=Thresholds,",(job-1)*b_size,",DBpath='",DBpath,"',name='",name,"')\n",sep="")),file=qtlfile,append=T)
	cat("report(3,\"JobFinished\")\n",file=qtlfile,append=T)
	#Quit
	cat("q(\"no\")","\n",sep="",file=qtlfile,append=T)
}

generateESTfile <- function(DBpath = "",nperm,map,method,model,taskID){
	#Generates a qtlfile to estimate walltime
	qtlfile <- paste("~/run",taskID,"/ESTtime.R",sep="")
	#Print our report function
	cat("\nreport <- function(status,text){\n",file=qtlfile)
	cat("\ttask <- ",taskID,"\n",file=qtlfile,append=T)
	cat("\ttext <- substr(URLencode(text),0,100)\n",file=qtlfile,append=T)
	cat("\tlink <- paste(\"",DBpath,"/taskreporter?job=\",task,\"&subjob=0&statuscode=\",status,\"&statustext=\",text,sep=\"\")\n",sep="",file=qtlfile,append=T)
	cat("\tgetURL(link)\n",file=qtlfile,append=T)
	cat("\tif(status==-1){\n\t\tcat(\"!!!\",text,\"!!!\")\n\t\t\n\t\tq(\"no\")\n\t}\n",file=qtlfile,append=T)
	cat("}\n\n",file=qtlfile,append=T)
	#load needed libraries	
	cat("library(qtl,lib.loc='~/libs')","\n",sep="",file=qtlfile,append=T)
	cat("library(bitops,lib.loc='~/libs')","\n",sep="",file=qtlfile,append=T)
	cat("library(RCurl,lib.loc='~/libs')","\n",sep="",file=qtlfile,append=T)
	cat("source(\"",paste(dbpath,"/api/R",sep=""),"\")\n",file=qtlfile,append=T)
	cat(Generate_Statement(paste("load(\"~/run",taskID,"/cross.RData\")","\n",sep="")),file=qtlfile,append=T)
	cat("cross$pheno <- cross$pheno[1]","\n",sep="",file=qtlfile,append=T)
	cat(Generate_Statement(paste("results <- bootstrap(cross,n.run=",nperm,",Funktie=",map,",method='",method,"',model='",model,"',multicore=FALSE,verbose=TRUE)","\n",sep="")),file=qtlfile,append=T)
	#Quit
	cat("q(\"no\")","\n",sep="",file=qtlfile,append=T)
}
