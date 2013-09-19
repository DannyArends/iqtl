#
#
# geneontology.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified jun, 2011
# first written jun, 2011
# 
# Gene ontology routines for eQCL Analysis
# geneMergeGO using the geneMerge perl script
#

geneMergeGO <- function(input, annotation, filename = "output", peekheight=1, inputThreshold=0.25, type=c("CC","BP","MF"), geneMerge = "E:/GeneMerge1.2-Windows/GeneMerge1.2.pl", logfile){
  abovecutoff <- which(input[[4]] >= peekheight)
  cat("Found: ",length(abovecutoff), "peeks at:", peekheight, ",", inputThreshold, "\n", file=logfile, append=TRUE)
  while(length(abovecutoff) >= 15){
    peekheight = peekheight+1
    abovecutoff <- which(input[[4]] > peekheight)
  }
  resFile <- sub("csv", "txt", paste("analysis/peek_summary_",filename,sep=""))
  cat("Found: ", length(abovecutoff), "peeks at:", peekheight, ",", inputThreshold, "\n", file=logfile, append=TRUE)
  if(length(abovecutoff) > 0){
    for(seed in names(abovecutoff)){
      cat("Analyzing:", seed, "\n")
      cat("Analyzing:", seed, "\n", file=logfile, append=TRUE)
      underpeek <- names(which(abs(input[[1]][seed,]) > inputThreshold))
      if(length(underpeek) > 10){
        cat(annotation[underpeek],sep="\n",file=paste("analysis/tmp",type,".txt",sep=""))
        system(paste("perl ",geneMerge," yeast.",type,"  yeast.",type,".use Uni2500.txt analysis/tmp",type,".txt analysis/out",type,".txt",sep=""),FALSE,show.output.on.console=FALSE)
        GOdata <- read.csv(paste("analysis/out",type,".txt",sep=""),sep="\t",colClasses="character")
        cat("GO done for ",seed,"-",underpeek,"\n",file=logfile,append=TRUE)
        if(dim(GOdata)[1] !=0){
          cat("GO significant!!!\n",file=logfile,append=TRUE)
          for(x in 1:dim(GOdata)[1]){
            cat(seed, annotation[seed], attr(input,"marker"), type, as.character(GOdata[x,]),"\n",file=resFile,sep="\t",append=TRUE)
          }
        }
        file.remove(paste("analysis/tmp",type,".txt",sep=""))
        file.remove(paste("analysis/out",type,".txt",sep=""))
      }else{
        cat("Skipped: ",seed,"eQCL peek to low:",length(underpeek),"\n",file=logfile,append=TRUE)
      }
      gc()
    }
  }
}
