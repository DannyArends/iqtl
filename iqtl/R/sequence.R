#
#
# sequence.R
#
# copyright (c) 2011 Danny Arends
# last modified Mrt, 2011
# first written Mrt, 2010
# 
# R functions to do load sequences from fasta format
#

#reads a single fasta file specified by directory and filename
read.fasta <- function(directory, filename){
  cat(paste(directory,"/",filename,sep=""),"\n")
  rawsequence <- read.table(paste(directory,"/",filename,sep=""),skip=1,sep="")
  res <- NULL
  for(x in 1:nrow(rawsequence)){
    res <- paste(res,as.character(rawsequence[x,]),collapse="",sep="")
  }
  res
}

#reads in all the files in the directory
#returns a list with at each element the sequence of the chromosome
read.genome <- function(directory="genome_sequence"){
  files <- dir(directory)
  genome <- vector("list",length(files))
  cnt <- 1
  for(x in files){
    genome[[cnt]] <- read.fasta(directory,x)
    attr(genome[[cnt]],"name") <- x
    cnt <- cnt + 1
  }
  genome
}
