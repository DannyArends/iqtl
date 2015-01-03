######################################################################
# snow.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: mqmmultitoscanone
#
#
######################################################################

create.client <- function(clustersetup=NULL, host=NULL, rscript=NULL, snowlib=NULL, cores=1, verbose=TRUE){
  if(is.null(clustersetup)){
    if(verbose) cat("Starting new clustersetup\n")
  }else{
    if(verbose) cat("Adding new machine to cluster\n")
  }
  client <-list(host = host, rscript = rscript, snowlib = snowlib)
  clustersetup <- c(clustersetup,rep(client,cores))
  clustersetup  
}

create.cluster <- function(clustersetup=NULL,verbose=TRUE){
  if(!is.null(clustersetup)){
    if(verbose) cat("Starting new cluster\n")
    cl <- makeCluster(clustersetup,type = "SOCK",port=3456,outfile="d:/test.txt",master="localhost")
		if(verbose) cat("Starting new clustercall\n")
    clusterCall(cl, function() Sys.info()[c("nodename","machine")])
    if(verbose) cat("Stopping cluster\n")
    stopCluster(cl)
  }
}


client.start <- function(verbose=TRUE){
  connection <- socketConnection(Sys.info()["nodename"], port = 3456,server=T,blocking=T)
  monitoring<-TRUE
  
  targetfilename <- tempfile()
  targetfile <- file(targetfilename)
  while(isOpen(connection) && monitoring) {
    Sys.sleep(1); 
    command <- readLines(connection,n=1)
    commandforclient<-FALSE
    if(command=="#QUIT"){
      if(verbose) cat("Received the command to go down\n")
      close(targetfile)
      monitoring=FALSE
      commandforclient=TRUE
    }
    if(command=="#NEWFILE"){
      if(verbose) cat("Received the command to build a new file\n")
      close(targetfile)
      targetfilename <- tempfile()
      targetfile <- file(targetfilename)
      commandforclient=TRUE
    }
    if(command=="#EXECUTE"){
      if(verbose) cat("Received the command to execute the current file\n")
      cat(file=targetfile,"q('no')\n",append=T)
      close(targetfile)
      setwd(tempdir())
      system(paste("R CMD BATCH",targetfilename," --vanilla"))
      targetfilename <- tempfile()
      targetfile <- file(targetfilename)
      commandforclient=TRUE
    }
    if(!commandforclient) cat(file=targetfile,command,"\n",append=T)
    if(verbose) cat(command,"\n")
  }
  close(connection);
}

test.cluster <- function(){
  macrscript.default = "/Library/Frameworks/R.framework/Resources/bin/Rscript"
  snowlibmac.default = "/Library/Frameworks/R.framework/Resources/library/snow"

  linuxrscript.default = "/usr/lib64/R/bin/Rscript"
  linuxsnowlib.default = "/home/luke/tmp/lib"
  
  winrscript.default = "C:/Program Files/R/R-2.11.0/bin/Rscript.exe"
  winsnowlib.default = "C:/Program Files/R/R-2.11.0/library/"

  clustersetup <- create.client(host="localhost",rscript=winrscript.default,snowlib=winsnowlib.default)
  #clustersetup <- create.client(host="129.125.143.25",rscript=winrscript.default,snowlib=winsnowlib.default,clustersetup=clustersetup)
  create.cluster(clustersetup)
  client.start()
}
