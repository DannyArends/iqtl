#####################################################################
#
# Realbasic.R
#
# copyright parent_check (c) 2009, Danny Arends
# last modified Apr, 2009
# first written Apr, 2009
# 	ourcat	CAT with verbose switch
#	Rn		Randomnumber
#	RnV		Vector of Randomnumbers
#	RnM		Matrix of Randomnumber
#	Names	Vector of names with PrefixNum
#
######################################################################

ourcat <- function(file=FALSE,name="TEST",...,verbose=TRUE,append=TRUE){
	if(verbose){
		cat(...)
	}
	if(file){
		if(name == "RESULT"){
			filename <- paste("Results_",format(Sys.time(), "%d%b%Y_%Hh%M"),".txt",sep="")
		}else{
			filename <- paste(name,format(Sys.time(), "%d%b%Y_%Hh%M"),"log.txt",sep="")
		}
		cat(file=filename,...,append=append)
	}
}

Rn <- function(num){
	as.integer(runif(1)*(num+1))
}

RnN <- function(num){
	sum(replicate(num,Rn(1)))
}

RnV <- function(num,times){
	replicate(times,Rn(num))
}

RnVN <- function(num,times){
	replicate(times,RnN(num))
}

RnM <- function(num,col,row){
	replicate(col,rbind(RnV(num,row)))
}

RnMN <- function(num,col,row){
	replicate(col,rbind(RnVN(num,row)))
}

myNames <- function(prefix,num){
	paste(prefix,rep(1:num),sep="")
}

smoothV <- function(data,smoothF,num){
	result <- NULL
	for(x in 1:length(data)){
		begin <- x-(num/2)
		end <- x+(num/2)
		if(begin > 1 && end < length(data)){
		result <- c(result,smoothF(as.numeric(data[begin:end])))
		}
		if(begin < 1){
		result <- c(result,smoothF(as.numeric(data[1:end])))
		}
		if(end > length(data)){
		result <- c(result,smoothF(as.numeric(data[begin:length(data)])))
		}
	}
	result
}

Names <- function(string="I",number=1){
  res <- paste(string,1:number,sep="")
  res
}
