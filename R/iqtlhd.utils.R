#
# iqtlhd.utils.R
#
# copyright (c) 2012 Danny Arends
# last modified Jan, 2012
# first written Jan, 2012
# 

connected <- function(){
  output <- (.C("connected",x=as.integer(0), PACKAGE="qtlHD")$x==1)
  cat("qtlHD:",output,"\n")
  invisible(output)
}
