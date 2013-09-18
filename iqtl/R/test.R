

#library(iqtl)
#qtl <- read.table("qtl_2500.txt")
#difcount <- read.table("difcountmatrix_2500.txt")
#cross <- read.cross("csvr",file="yeast_2500.csv",geno=c(0,1))

#cnt <- 0
#for(x in 1:ncol(qtl)){
#  if(max(difcount[,x]) > 15){
#    q <- lodscorestoscanone(cross,qtl[,x])
#    d <- lodscorestoscanone(cross,difcount[,x])
#    jpeg(paste("analysis/",colnames(qtl)[x],".jpg",sep=""))
#    plot(q,d)
#    dev.off()
#    cnt <- cnt+1
#  }
#}
