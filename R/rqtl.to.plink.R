#
# rqtl.to.plink.R
#
# copyright (c) 2011 Danny Arends and Ritsert C. Jansen
# last modified Jan, 2012
# first written Jan, 2012
#
# get.chr, rqtl.to.plink
#

get.chr <- function(cross,m=5){
  chr <- 1
  nmar <- length(pull.map(cross)[[chr]])
  while(m > nmar){
    m <- m - nmar
    chr <- chr + 1
    nmar <- length(pull.map(cross)[[chr]])
  }
  chr;
}

rqtl.to.plink <- function(cross, phenotype=1){
  #.PED
  if("sex" %in% colnames(hyper$pheno)){
    sex <- cross$pheno$sex
    sex <- sex=="male"
    sex[which(sex==TRUE)] <- 1
    sex[which(sex==FALSE)] <- 2
  }else{
    sex <-  rep(9,nind(cross))
  }
  genotypes <- pull.geno(cross)
  if(any(genotypes==0,na.rm=T)) genotypes <- genotype+1
  genotypes[is.na(genotypes)] <- 0 #0 is NA
  for(ind in 1 : nind(cross)){
    cat(1,ind, 0, 0, sex[ind], cross$pheno[ind, phenotype],as.numeric(genotypes[ind,]),"\n", file="cross.ped")
  }
  cat("Wrote PED file",nind(cross),"individuals at",sum(nmar(cross)),"loci\n")
  #.MAP
  for(m in 1:sum(nmar(cross))){
    cat(get.chr(cross,m), markernames(cross)[m], unlist(pull.map(cross))[m], unlist(pull.map(cross))[m]*1000000,"\n", file="cross.map")
  }
  cat("Wrote MAP file at",sum(nmar(cross)),"loci\n")
}

# end of rqtl.to.plink.R
