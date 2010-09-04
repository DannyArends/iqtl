#####################################################################
#
# common.R
#
# copyright (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Nov, 2009
# 
# Common functions of the ClusterJobs package
# Contains: run_cluster
#
######################################################################

######################################################################
#
# generateRunfile: Generates a .sh file with job description to be submitted to pbs
# generateQTLfile: Generates a QTL file which is executed (R CMD BATCH <qtlfile>) by a cluster
# DownloadnSave: Saves a Crossobject to the HDD in .RData format for easy loading
# run_cluster: main routine to distribute jobs across a cluster of computers
# generateRunfile: Generates a .sh file with job description to be submitted to pbs
# generateQTLfile: Generates a QTL file which is executed (R CMD BATCH <qtlfile>) by a cluster
# generateESTfile: Generates a QTL file to estimate runtime for each item
# mqmmultitomatrix: Generates a matrix, that can be uploaded by using the add.datamatrix function
#
######################################################################

mqmmultitomatrix <- function(mqmmulti){
  matr <- NULL
  coln <- NULL
  for(x in 1:length(mqmmulti)){
    matr <- cbind(matr,mqmmulti[[x]][,3])
    coln <- c(coln, substr(colnames(mqmmulti[[x]])[3],5,nchar(colnames(mqmmulti[[x]])[3])))
  }
  rownames(matr) <- rownames(mqmmulti[[1]])
  colnames(matr) <- coln
  matr
}

generateRunfile <- function(job, est,jobid){
	#Generate a runfile to submit to the cluster
	runfile <- paste("~/run",jobid,"/run",job,".sh",sep="")
	cat("#!/bin/sh","\n",sep="",file=runfile)
	#We need just 1 node TODO: use 2 processors and SNOW on a clusterNODE
	cat("#PBS -l nodes=1","\n",sep="",file=runfile,append=T)
	#Our estimate
	cat("#PBS -l walltime=",est,"\n",sep="",file=runfile,append=T)
	#Go to the coreect location
	cat(paste("cd $HOME/run",jobid,sep=""),"\n",sep="",file=runfile,append=T)
	#StartRunning
	cat("R CMD BATCH subjob",job,".R",sep="",file=runfile,append=T)
	runfile
}

DownloadnSave <- function(investigationname,DBmarkerID = "", DBtraitID = "", dbpath = "",jobid,njobs){
	#Generates a R-script to download all the information and build a cross object
	qtlfile <- paste("~/run",jobid,"/download.R",sep="")
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
    cat("source(\"",paste(dbpath,"/api/R",sep=""),"\")\n",sep="",file=qtlfile,append=T)
	#Downloading of Cross object (secured)
	cat(Generate_Statement(paste("cross <- CrossFromMolgenis(genotypematrixname='",DBmarkerID,"',phenotypematrixname='",DBtraitID,"',investigationname='",investigationname,"')","\n",sep="")),file=qtlfile,append=T)
	cat(Generate_Statement(paste("save(cross,file=\"~/run",jobid,"/cross.RData\")","\n",sep="")),file=qtlfile,append=T)
	cat("q(\"no\")","\n",sep="",file=qtlfile,append=T)
}


run_cluster <- function(name = "qtlTEST",investigation="",genotypes = "", phenotypes = "",totalitems=1, njobs=1,dbpath = "",jobid=1,map="scanall",method="hk",model="normal",step=0){
	#Initializes the cluster, installs R-libraries (could only be done once)
	#Downloads the datafile from molgenis and save it as a .RData
	#Estimated time needed for each run (very crude) (totalitems in run + 25%)
	#--Generate a QTLfile
	#--Generate runfile for cluster
	#--Sends the runfile as a job to the cluster
	totalitems <- as.integer(totalitems)
	njobs <- as.integer(njobs)
	jobid <- as.integer(jobid)	
	prepare_cluster(jobid)
	cat("000000000000000000000000000000\n")
	library(qtl,lib.loc='~/libs')
	cat("1\n")
	library(bitops,lib.loc='~/libs')
	library(RCurl,lib.loc='~/libs')
	cat("2\n")
	source(paste(dbpath,"/api/R",sep=""))
	cat("3\n")
	est = NULL
	runfile = NULL
	if(totalitems < njobs){
		njobs = totalitems
	}
	nprun <- ceiling(totalitems/njobs)
	cat("# of Traits:",totalitems,"\n# of Jobs",njobs,"# per run",nprun,"\n")
	tryCatch(DownloadnSave(investigationname=investigation, genotypes,phenotypes,dbpath=dbpath,jobid,njobs)
		,error =  function(e){report(dbpath,jobid,0,-1,"Downloadscript")}
	)
	report(dbpath,jobid,0,2,"GeneratedDownload")
	tryCatch(system("R CMD BATCH download.R")
		,error =  function(e){report(dbpath,jobid,0,-1,"DownloadingCrossobject")}
	)
	report(dbpath,jobid,0,2,"FinishedDownloadingDatasets")	
	tryCatch(est <- est_runtime(njobs,totalitems,dbpath=dbpath,nprun,jobid,map,method,model,step)
		,error =  function(e){report(dbpath,jobid,0,-1,"EstimatingTime")}
	)
	if(is.null(name)) name="qtlTEST"
	cat("Estimated for ",name,". Runtime per job=",est,"\n")
	report(dbpath,jobid,0,2,"EstimatedRuntime")
	#ALL DONE NOW WE CAN GO INTO/run directory and make some calculations
	for(x in 1:njobs){
		cat("Generating: ",x,".1/",njobs,"\n",sep="")
		tryCatch(generateQTLfile(dbpath=dbpath, job=x, b_size=nprun,Ntraits=totalitems,name=name,jobid,map,method,model,step,investigationname=investigation)
			,error =  function(e){report(dbpath,jobid,x,-1,"GeneratingQTLrunfile")}
		)
		report(dbpath,jobid,0,2,paste("Generated_QTL",x,sep=""))	
		cat("Generating: ",x,".2/",njobs,"\n",sep="")
		tryCatch(runfile <- generateRunfile(x,est,jobid)
			,error =  function(e){report(dbpath,jobid,x,-1,"GeneratingSHrunfile")}
		)
		report(dbpath,jobid,0,2,paste("Generated_SH",x,sep=""))	
		cat("Submitting: ",x,".3/",njobs,":",runfile,"\n",sep="")
		#OLD call to sh to compute on the Sceduler
		tryCatch(system(paste("qsub ",runfile,sep=""))
			,error =  function(e){report(dbpath,jobid,x,-1,"SubmissionPBS:")}
		)
		report(dbpath,jobid,0,2,paste("Submitted_",x,sep=""))	
		report(dbpath,jobid,x,1,"TaskqueuedbyPBS")
		#system(paste("sh ",runfile,sep=""))
	}
	report(dbpath,jobid,0,3,"PreparationDone")
}

Generate_Statement <- function(statement){
	#printing a secured statement (just means we report a 3 if we fail)
	secured <- paste("tryCatch(\n\t",statement,"\t,error =  function(e){report(-1,e$message)}\n)\n",sep="")
	secured
}

est_runtime <- function(njobs,ntraits=1, dbpath = "",num_per_run=1,jobid,nperm,map,method,model){
	#Crude estimation of time that it would take a job of num_per_run qtls to finish, we get all the data from molgenis and run 2 qtls profiles
	generateESTfile(dbpath,nperm,map,method,model,jobid)
	#time execution of executing 1 trait
	s <- proc.time()
	system("R CMD BATCH ESTtime.R")
	e <- proc.time()
	#Add 5% for security
	EST <- (num_per_run)*((e[3]-s[3])+10) 
	ESTtime <- sprintf("%02.f:%02.f:%02.f",EST %/% 3600, (EST%%3600) %/% 60, round(EST%%60, digits = 0))
	if(EST >  864000){
		report(dbpath,jobid,0,-1,"JobsTooLong")
	}else{
		ESTtime
	}
}

install_libraries <- function(){
	cat("Installing libraries\n")
	#Commands to execute to install R-libs
	system("R CMD INSTALL ~/required/qtl_1.15-1.tar.gz --library=~/libs")
	system("R CMD INSTALL ~/required/RCurl_0.91-0.tar.gz --library=~/libs")
	system("R CMD INSTALL ~/required/snow_0.3-3.tar.gz --library=~/libs")
	system("R CMD INSTALL ~/required/ClusterJobs_0.1.tar.gz --library=~/libs")
}

prepare_cluster <- function(jobid){
	cat("Creating directory\n")
	#Create a directory and switch to it (in R and on the shell (so our next R's we'll spawn will not have 2 switch dirs)
	system(paste("mkdir run",jobid,sep=""))
	system(paste("cd run",jobid,sep=""))
	setwd(paste("~/run",jobid,sep=""))
}

report <- function(dbpath,task,job,status,text){
	progress <- 0
	text <- substr(URLencode(text),0,100)
	link <- paste(dbpath,"/taskreporter?job=",task,"&subjob=",job,"&statuscode=",status,"&statustext=",text,"&statusprogress=",progress,sep="")
	getURL(link)
	if(status==-1){
		q("no")
	}
}
