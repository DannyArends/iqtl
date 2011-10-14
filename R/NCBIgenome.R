#Basic search NCBI function
searchNCBIgenome <- function(species = "", verbose = TRUE){
  require("RCurl")
  ftplocation <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/"
  querystring <- paste(ftplocation,species,"/",sep="")
  if(verbose) cat(querystring,"\n",sep="")
  HTMLdata <- getURL(querystring,ftp.use.epsv = FALSE, dirlistonly = TRUE)
  formatted <- parseFTPdir(HTMLdata)
  if(species != ""){
    chr <- formatted[which(substr(formatted,0,3)=="CHR"),]
    if(length(chr) > 0){
      CHRdata <- vector("list",length(chr))
      CHRcnt <- 1
      for(x in chr){
        CHRdata[[CHRcnt]] <- paste(ftplocation,species,"/",x,"/",searchNCBIgenome(paste(species,"/",x,sep="")),sep="")
        CHRcnt <- CHRcnt + 1
      }
      names(CHRdata) <- chr
      CHRdata
    }else{
      formatted
    }
  }else{
    formatted
  }
}

#Parse the format of an FTP directory listing
parseFTPdir <- function(FTPdata){
 regexp_p <- "(.+?)\r\n"
 possible <- regexpr(regexp_p, FTPdata)
 mydata <- NULL
  while(possible > 0){
    mydata <- rbind(mydata,substr(FTPdata,possible[[1]],possible[[1]]+(attr(possible,"match.length")-3)))
    FTPdata <- substr(FTPdata,possible[[1]]+attr(possible,"match.length"),nchar(FTPdata))
    possible <- regexpr(regexp_p, FTPdata)
  }
  invisible(mydata)
}

#Download the files to disk
getNCBIfiles <- function(species, extension=".faa", verbose = TRUE){
  for(x in 1:length(species)){
    toDownload <- which(extension == substr(species[[x]],nchar(species[[x]])-3,nchar(species[[x]])))
    if(verbose) cat("Downloading:",species[[x]][toDownload],"\n")
    filedata <- getURL(species[[x]][toDownload])
    if(verbose) cat("Saving:",paste(getwd(),names(species)[x],extension,sep=""),"\n")
    cat(filedata,file=paste(names(species)[x],extension,sep=""))
  }
}

#Checks if your species is available at NCBI genome
hasSpecies <- function(name = "Arabidops", species, ignore.case = TRUE){
  possible <- grep(name, species, ignore.case=ignore.case)
  if(length(possible) > 1){
    warning(paste("Multiple found, did you mean: ",paste(species[possible],collapse=", or "),"\n",sep=""))
    return(species[possible])
  }
  if(length(possible) == 1){
    return(species[possible])
  }
  stop(paste("No species called '",name,"' is not found",sep=""))
}

test.NCBIgenome <- function(){
  species <- searchNCBIgenome()
  #Multiple possible matches, issues a warning
  hasSpecies("arab",species)
  #Model plants/animals
  AT <- searchNCBIgenome(hasSpecies("thaliana",species))
  ELEGANS <- searchNCBIgenome(hasSpecies("elegans",species))
  #Human datasets
  HUMAN <- searchNCBIgenome(hasSpecies("elegans",species))
  getNCBIfiles(AT)
}
