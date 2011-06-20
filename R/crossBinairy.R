
getMarkerAnnotation <- function(cross){
  distance <- as.numeric(unlist(pull.map(cross)))
  chr <- rep(names(cross$geno),nmar(cross))
  markername <-  as.character(unlist(lapply(FUN=names,pull.map(cross))))
  markertype <- as.character(rep(unlist(lapply(cross$geno,function(x){attr(x,"class")})),nmar(cross)))
  data.frame(markername,chr,distance,markertype)
}


test.binary <- function(cross,name="test.bin"){
  genotypes <- pull.geno(cross)
  phenotypes <- pull.pheno(cross)
  markerannotation <- getMarkerAnnotation(cross)
  
  #store first one
  .C("r_writebinaryfile",
    name=as.character(name),            #file name
    nlength=as.integer(nchar(name)),    #namelength
    nrows=as.integer(nind(cross)),      #rowdimension matrix
    ncols=as.integer(sum(nmar(cross))), #coldimension matrix
    type=as.integer(1),                 #type of data
    lengths=as.integer(1),              #lengths of dataelements
    mdata=as.character(genotypes)       #data
  )
}

write.cross.binary <- function(cross,name="test.bin"){
  genotypes <- pull.geno(cross)
  phenotypes <- pull.pheno(cross)
  markerannotation <- getMarkerAnnotation(cross)
  
  #store first one
  .C("r_writebinaryfile",
    name=as.character(name),            #file name
    nlength=as.integer(nchar(name)),    #namelength
    nrows=as.integer(nind(cross)),      #rowdimension matrix
    ncols=as.integer(sum(nmar(cross))), #coldimension matrix
    type=as.integer(3),                 #type of data
    lengths=as.integer(1),              #lengths of dataelements
    mdata=as.character(genotypes)       #data
  )

  #store second one  
  .C("r_writebinaryfile",
    name=as.character(name),          #file name
    nlength=as.integer(nchar(name)),  #namelength
    nrows=as.integer(nind(cross)),    #rowdimension matrix
    ncols=as.integer(nphe(cross)),    #coldimension matrix
    type=as.integer(2),               #type of data
    lengths=as.integer(1),            #lengths of dataelements (IGNORED HERE)
    mdata=as.double(phenotypes)       #data
  )
  
  #store third one  
  .C("r_writebinaryfile",
    name=as.character(name),              #file name
    nlength=as.integer(nchar(name)),      #namelength
    nrows=as.integer(nind(cross)),        #rowdimension matrix
    ncols=as.integer(4),                  #coldimension matrix
    type=as.integer(4),                   #type of data
    lengths=as.integer(1),                #lengths of dataelements
    mdata=as.character(markerannotation)  #data
  )
}