#
#
# gnuplot.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: creategnuplot, createdatafile, testgnuplot, createprofilebymarkerscript, createheatmapbymarkerscript
#
#

creategnuplot <- function(cross,result,output=c("screen","svg","postscript"),type=c("heatmap","surface","lines","points","linespoints","dots","impulses"),pheno.col=1:length(result)){
  Tfile <- file("plot.p", "w")
  cat("#Scriptfile\n",file=Tfile)
  if(output[1]=="svg"){
    cat("set terminal svg\n",file=Tfile)
    cat("set output \"outputfile.svg\"\n",file=Tfile)
  } 
  if(output[1]=="postscript"){
    cat("set terminal postscript color solid\n",file=Tfile)
    cat("set output \"outputfile.ps\"\n",file=Tfile)
  }  
  for(x in pheno.col){
    createdatafile(result,type="2D",pheno.col=x)
  }
  createdatafile(result,type="3D",pheno.col=pheno.col)
  if(type[1]!="heatmap" && type[1]!="surface"){
    cat("#Phenotype ",pheno.col,"\n",sep=",",file=Tfile)
    createprofilebymarkerscript(result,type=type,pheno.col=pheno.col,Tfile=Tfile)
  }else{
    cat("#Phenotypes ",pheno.col,"\n",sep="",file=Tfile)  
    cat("set autoscale\nunset log\nunset key\nunset label\n",file=Tfile)
    createheatmapbymarkerscript(result,type=type,pheno.col=pheno.col,Tfile=Tfile);
  }
  cat("pause -1\n",file=Tfile)
  close(Tfile)
}


createprofilebymarkerscript <- function(result,type=c("lines","points","linespoints","dots","impulses"),pheno.col,Tfile,size=c(1.0,1.0),loc=c(0.0,0.0),grid=FALSE){
  temp <- lapply(result,getThird)
  temp <- do.call("rbind",temp)
  
  #cat("set size ",size[1],",",size[2],"\n",file=Tfile)
  #cat("set origin ",loc[1],",",size[2],"\n",file=Tfile)
  if(grid) cat("set grid\n",file=Tfile)
  cat("set title \"Profile - Scale by Marker\"\n",file=Tfile)
  cat("set xlabel \"Marker\"\n",file=Tfile)
  for(x in 1:nrow(result[[1]])){
    cat("set label ",x," \"",rownames(result[[1]])[x],"\" at ",x,",0 rotate by -90 left\n",sep="",file=Tfile)
  }
  cat("set ylabel \"LOD\"\n",file=Tfile)
  cat("set palette rgbformula 7,-8,-3\n",file=Tfile)
  cat("set yrange [",round(min(temp),digits=0),":",round(max(temp),digits=0),"]\n",sep="",file=Tfile)
  cat("set xrange [0.5:",nrow(result[[1]]),".5]\n",sep="",file=Tfile)
  cat("set style fill transparent pattern 4 bo\n",sep="",file=Tfile)
  cat("plot ",sep="",file=Tfile)
  count <- 1
  for(x in pheno.col){
    if(count >1) cat(", ",sep="",file=Tfile)
    cat("\"data",x,".dat\" with ",type," title \"",names(pull.pheno(result))[x],"\" ",sep="",file=Tfile)
    count <- count +1
  }
   cat("\n",sep="",file=Tfile)
}

createheatmapbymarkerscript <- function(result,type=c("heatmap","contour"),pheno.col,Tfile,size=c(1.0,1.0),loc=c(0.0,0.0),grid=FALSE){
  temp <- lapply(result,getThird)
  temp <- do.call("rbind",temp)
  
  #cat("set size ",size[1],",",size[2],"\n",file=Tfile)
  #cat("set origin ",loc[1],",",size[2],"\n",file=Tfile)
  if(grid) cat("set grid\n",file=Tfile)

  if(type[1]=="heatmap"){
    cat("set title \"Heatmap - Scale by Marker\"\n",file=Tfile)
    cat("set xlabel \"Marker\"\n",file=Tfile)
    cat("set ylabel \"Trait\"\n",file=Tfile)
    cat("set ytics (",paste(paste("\"",names(pull.pheno(result))[pheno.col],"\"",collapese="")," ",seq(0.5,24.5,step=1),collapse=","),")\n",file=Tfile,sep="")
    cat("set palette rgbformula 7,-8,-3\n",file=Tfile)
    cat("set xrange [0.5:",nrow(result[[1]]),".5]\n",sep="",file=Tfile)
    cat("set yrange [0.5:",length(result),".5]\n",sep="",file=Tfile)  
    cat("set cbrange [",round(min(temp),digits=0),":",round(max(temp),digits=0),"]\n",sep="",file=Tfile)
    cat("set cblabel \"LOD Score\"\n",file=Tfile)
    cat("set view map\n",file=Tfile)
    cat("splot \"data.dat\" using 2:1:3 with image\n",file=Tfile)
  }
  if(type[1]=="surface"){
    cat("set title \"Surface - Scale by Marker\"\n",file=Tfile)
    cat("set ylabel \"Marker\"\n",file=Tfile)
    cat("set xlabel \"Trait\"\n",file=Tfile)
    cat("set xtics (",paste(paste("\"",names(pull.pheno(result))[pheno.col],"\"",collapese="")," ",seq(0.5,24.5,step=1),collapse=","),")\n",file=Tfile,sep="")
    cat("set autoscale\n",file=Tfile)  
    cat("set pm3d at s\n",file=Tfile)
    cat("set palette gray\n",file=Tfile)
    cat("unset surface\n",file=Tfile)    
    cat("set zrange [",round(min(temp),digits=0),":",round(max(temp),digits=0),"]\n",sep="",file=Tfile)
    cat("set zlabel \"LOD Score\"\n",file=Tfile)
    cat("set parametric\n",file=Tfile)
    cat("set style data lines \n",file=Tfile)
    cat("splot \"data.dat\"\n",file=Tfile)
  }
}

createdatafile <- function(result,type=c("3D","2D"),pheno.col=1:length(result),by=c("cM","Marker")){
  if(type=="2D"){
    Tfile <- file(paste("data",pheno.col,".dat",sep="",collapse=""), "w")
  }else{
    Tfile <- file("data.dat", "w")
  }
  cat("#Datafile\n",file=Tfile)
  for(t in 1:length(result)){
    if(t %in% pheno.col){
      for(m in 1:nrow(result[[t]])){
        if(type[1]=="2D"){
          cat(c(m,round(result[[t]][m,3], digits = 1),"\n"),file=Tfile,append=T)
        }else{
          cat(c(t,m,round(result[[t]][m,3], digits = 1),"\n"),file=Tfile,append=T)
        }
      }
      cat("\n",file=Tfile,append=T)
    }
  }
  close(Tfile)
}

testgnuplot <- function(location="d:/generated",output=c("screen","svg","postscript"),...){
  if(!exists("result")){
    location <-"d:/generated"
    setwd(location)
    library(qtl)
    data(multitrait)                                              
    multitrait <- fill.geno(multitrait)
    result <- scanall(multitrait,scanfunction=mqmscan,logtransform=TRUE)  
  }
  creategnuplot(multitrait,result,output=output,...)
  system("gnuplot plot.p",wait=FALSE,invisible=FALSE,show.output.on.console=FALSE)
  result
}

#result <- testgnuplot()
#result <- testgnuplot(type="surface")
