#####################################################################
#
# crash_cluster.R
#
# copyright (c) 2009, Danny Arends
# last modified Sep, 2011
# first written Dec, 2009
# 
# Contains: endlessloop and support functions
#
# Designed to test clustersoftware for "crash-safty" by letting it 
# execute Rscripts (a commonly used language in bioinformatics) aimed 
# at bringing the executing OS down and thus destroying the running 
# GRID. Types of attacks:
#
# BUT FIRST: uncomment your type of attack:  ;) from mild to more severe
#
# type = "loop"             #just loop forever             Test if it can be detected
# type = "nocpu"            #loop no-cpu                   To fool programs just looking for cpu usage 
# type = "highcpu"          #loop Burning the cpu          Tests
# type = "mem"              #loop with increased memusage  Try n explode the surrounding shell, see if the grid can still stop it
# type = "threads"          #Spawning new threads          TODO: Try n explode the surrounding shell, see if the grid can still stop it
# type = "system"           #Spawning new system threads   Try n explode the surrounding shell, see if the grid can still stop it
#
# Submit via an PBS .sh file containing: R CMD BATCH crash_cluster.R
# endlessloop(type)
######################################################################

endlessloop <- function(type){
  cnt <- 0
  list <- NULL
  while(TRUE){                                  #omg lol
    switch(type,
      loop = loop(cnt),
      nocpu = nocpu(10),
      highcpu = list <- c(list,highcpu(list)),
      mem = list <- c(list,mem(list)),
      threads = threads(FALSE,cnt),
      system = threads(TRUE,cnt)
    )
    cat("Loop ",cnt," Finished\n")
    cnt <- cnt+1
  }
  #gctorture(on = FALSE)                        #Why turn it off =) its an attack
  list
}

loop <- function(cnt){
  #Just run like any other lame ass job
  id <- round(runif(1)*100)
  cat("Loop ",cnt," Reporting Sir\n")
  if(id==42){
    string <- "42\n ... The meaning of Life, The universe and Everything ...\n"
    for(x in 1:nchar(string)){
      cat(substr(string, 1, x),"\n")
    }
  }
}

nocpu <- function(time){
  Sys.sleep(time)
}

highcpu <- function(list){
  gctorture(on = TRUE)
  newlist <- c(list,1:1024)           
  #Sys.sleep(1)
  newlist
}

mem <- function(list){
  newlist <- c(list,1:1024)
  #Sys.sleep(1)
  newlist
}

threads <- function(system, cnt){
  if(system){
    myfile <- paste("mysleeper",cnt,".R",sep="")
    cat("Sys.sleep(1);q(\"n\")",cnt,".R",sep="",file=myfile)
    system(paste("R CMD BATCH ",myfile,sep=""),wait=FALSE)
    #Sys.sleep(1)
  }else{
    #TODO add internal c threader
  }
}
