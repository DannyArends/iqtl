#
# scanonetomqmmulti.R
#
# copyright (c) 2010, Danny Arends
# last modified dec, 2010
# first written dec, 2010
# 
# R functions: scanonetomqmmulti
#


#Changes a qtl matrix (rows: traits, cols: markers) to a mqmmulti object
qtlmatrixtomqmmulti <- function(cross,qtlmatrix){
  rlist <- vector("list", nrow(qtlmatrix))
  for(x in 1:nrow(qtlmatrix)){
    rlist[[x]] <- lodscorestoscanone(cross,qtlmatrix[x,])
  }
  rlist
}

#Changes a scanone object (with multiple qtl profiles) to a mqmmulti object
scanonetomqmmulti <- function(cross,scanone){
  rlist <- vector("list", ncol(scanone)-3)
  for(x in 3:ncol(scanone)){
    rlist[[x-2]] <-lodscorestoscanone(cross,scanone[,x])
  }
  rlist
}
