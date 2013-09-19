######################################################################
#
# utils.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: mqmmultitoscanone, mqmmodelsasnetwork
#
#
######################################################################

mqmmultitoscanone <- function(result){
  if(!any(class(result)=="mqmmulti")){
    stop("No mqmmulti object. Please supply a valid mqmmulti object.") 
  }
  obj <- result[[1]][1:2]
  m <- do.call(cbind,lapply(FUN=getThird,result))
  obj <- cbind(as.matrix(obj),m)
  obj <- as.data.frame(obj)
  colnames(obj) <- c("chr","pos (cM)",rep("lod",ncol(m)))
  class(obj) <- c("scanone",class(obj))
  obj
}

mqmsupportint <- function(result,marker){
  if(!any(class(result)=="scanone")){
    stop("No scanone object. Please supply a valid scanone object.") 
  }
  num <- which(rownames(result)==marker)
  if(is.na(num&&1)) stop(paste("Cannot find marker named:",marker))
  result <- result[result[num,1]==result[,1],]
  num <- which(rownames(result)==marker)
  if(result[num,3] > 1.5){
    drop <- 1.5
  }else{
    drop <- 0.75 * result[num,3]
  }
  temp <- which(result[,3] < result[num,3]-drop)
  min <- temp[which(temp < num)[length(which(temp < num))]]
  max <- temp[which(temp > num)[1]]
  if(is.na(min&&1)){
    min <- 1
  }
  if(is.na(max&&1)){
    max <- nrow(result)
  }
  list(result[min,],result[max,])
}

#Retrieves the number of lines in a file
#Uses a buffer for reading so we don't explode out of memory
numLines <- function(filename, batchsize=2000, verbose=TRUE){
  if(missing(filename)) stop("Please supply a filename")
  file.in <- file(filename, 'rt')
  line <- readLines(file.in, n=1)
  s <- proc.time()
  cnt <- 1
  batch <- 0
  while (length(line)) {
    line<- readLines(file.in, n=1)
    cnt <- cnt +1
    if(cnt %% batchsize == 0){
      e <- proc.time()
      batch <- batch+1
      if(verbose) cat("Done with batch ",batch,"=",cnt,"lines",e[3]-s[3], "seconds\n")
    }
  }
  close(file.in)
  if(verbose) cat("# of lines: ",cnt,"\n")
  cnt
}

gcLoop <- function(verbose=FALSE){
	before <- gc()[2,3]
	bf <- before
	after <- gc()[2,3]
	while(before!=after){
		before <- after
		after <- gc()[2,3]
	}
	if(verbose) cat("Cleaned up memory from:",bf,"to:",after,"\n")
}
