#####################################################################
#
# plotroutines.R
#
# copyright parent_check (c) 2009, Danny Arends
# last modified Sep, 2010
# first written Mrt, 2009
# 
# Contains: several plot functions
#
######################################################################

plot.histkinship <- function(x,ind1,...){
	result1 <- x[[1]]
	result2 <- x[[2]]
	op <- par(mfrow = c(2,1))
	b <- seq(0,1,0.01)
	hist(result1[ind1,],breaks=b,main=paste("Kinshipscores individual:",ind1),xlab="Parental scores",...)
	for(x in seq(2,10,2)){
		lines(cbind(rep(mean(result1[ind1,],na.rm=T)+(x*sd(result1[ind1,],na.rm=T)),ncol(result1)),seq(0,(ncol(result1)-1))),col=palette()[x],lwd=1,lty=2)
		lines(cbind(rep(mean(result1[ind1,],na.rm=T)-(x*sd(result1[ind1,],na.rm=T)),ncol(result1)),seq(0,(ncol(result1)-1))),col=palette()[x],lwd=1,lty=2)
	}
	hist(result2[ind1,],breaks=b,main=paste("Kinshipscores individual:",ind1),xlab="Parental scores",...)
	for(x in seq(2,10,2)){
		lines(cbind(rep(mean(result2[ind1,],na.rm=T)+(x*sd(result2[ind1,],na.rm=T)),ncol(result2)),seq(0,(ncol(result2)-1))),col=palette()[x],lwd=1,lty=2)
		lines(cbind(rep(mean(result2[ind1,],na.rm=T)-(x*sd(result2[ind1,],na.rm=T)),ncol(result2)),seq(0,(ncol(result2)-1))),col=palette()[x],lwd=1,lty=2)
	}			
	op <- par(mfrow = c(1,1))
}

plot.parentplot <- function(x, start=1, num=(ncol(x)-start),...){
	x <- x[[1]]
	op <- par(mfrow = c(1,1))
	res <- x
	clrs <- rainbow(num+1)
	names <- NULL
	plot(res[start,],type="n",xlab="individual",ylab="kinship",...)
	for(x in start:(start+num)){
		names <- c(names,paste("Ind",x,sep=""))
		lines(res[x,],type="l",col=clrs[(x-start)+1])
	}
	legend("topright",legend=names,lty=1,lwd=1,col=clrs)
}

plot.parentimage <- function(x, start=1, num=(ncol(x)-start), cutoffsib=3,cutoffpar=5,...){
	x <- x[[1]]
	op <- par(mfrow = c(1,1))
	if(!(num <2)){
		x <- x[start:(start+num),start:(start+num)]
		meani <- mean(x,na.rm=T)
		sdi <- mean(sd(x,na.rm=T))
		image(x=start:(start+num),y=start:(start+num),z=x,xlab="Individual",ylab="Individual",col=c("white","lightgray","gray","darkgray","green","blue"),breaks=c(0,meani-15*sdi,meani,meani+sdi,meani+cutoffsib*sdi,meani+cutoffpar*sdi,(meani+100*sdi)),...)
		if(num <100){
			grid((num+1),lty="solid")
		}
		box()
	}else{
		stop("num needs to be larger than 2")
	}
}

plot.combined <- function(x,...){
	POPresult <- x[[1]]
	PARresult <- x[[2]]
	op <- par(mfrow = c(2,1))
	plot.popstats(POPresult,plots=c(F,T,F,T),...)
	plot.genotypes(PARresult)	
}

plot.all <- function(x,...){
	POPresult <- x[[1]]
	PARresult <- x[[2]]
	layout(matrix(c(1,1,3,2,2,2,4,4,4), 3,3, byrow = TRUE))
	op <- par(mfrow = c(2,1))
	plot.popstats(POPresult,plots=c(T,T,T,T),...)
	plot.genotypes(PARresult)	
}


plot.genotypes <- function(x,...){
	if(is.list(x)){
		op <- par(mfrow = c(1,1))
		x <- x[[2]]
	}
	ch3 <- c(0.0,1.0,7.82,11.35,16.27,25.0)
	colors <- gray(4:0/4)
	image(x=1:ncol(x),y=1:nrow(x),z=t(as.matrix(x)),axes=F,xlab="marker",ylab="Matingtype", col=colors, breaks=ch3,...)
	Axis(side=1,at=1:ncol(x),labels=1:ncol(x))
	Axis(side=2,at=1:nrow(x),labels=rownames(x))
}

plot.popstats <- function(x,smooth=F,num=10,smoothF=median,plots = c(T,T,T,F),...){
	if(is.list(x)){
		result <- x[[1]]
	}else{
		result <- x
	}
	ch1 <- c(3.84,6.64,10.83)
	PopulationScores <- result[,1]
	PopulationScores2 <- NULL
	PopulationScoressmooth <- NULL
	ParentalScores <- result[,4]
	ParentalScoressmooth <- NULL
	InbreadingScores <- result[,3]
	for(x in 1:nrow(result)){
		if(InbreadingScores[x] < 0){
			PopulationScores2 <- c(PopulationScores2,-as.numeric(PopulationScores[x]))
		}else{
			PopulationScores2 <- c(PopulationScores2,as.numeric(PopulationScores[x]))
		}			
	}
	if(smooth){
		PopulationScoressmooth <- smoothV(PopulationScores2,smoothF,num)
		ParentalScoressmooth <- smoothV(ParentalScores,smoothF,num)
	}
	op <- par(mfrow = c(sum(plots),1))
	#plot 1
	if(plots[1]){
		max <- max(as.numeric(PopulationScores2))
		min <- min(as.numeric(PopulationScores2))
		plot(c(1,length(PopulationScores)),c(min,max),xlab="markers",ylab="Chi^2",type='n',...)
		lines(PopulationScores2,lwd=2,col=rgb(0.3,0.3,1.0,1.0),type="o")
		if(smooth){
			lines(PopulationScoressmooth,lwd=2,col=rgb(0.0,0.0,0.0,1.0),type="l")
		}
		lines(rep(3.84,length(PopulationScores)),type="l",col=rgb(0.0,0.0,1.0,1.0),lty=2)
		lines(rep(-3.84,length(PopulationScores)),type="l",col=rgb(0.0,0.0,1.0,1.0),lty=2)
	}
	#plot 2
	if(plots[2]){
		max <- max(as.numeric(ParentalScores))
		min <- min(as.numeric(ParentalScores))
		plot(c(1,length(PopulationScores)),c(min,max),xlab="markers",ylab="Chi^2",type='n',...)
		lines(ParentalScores,lwd=2,col=rgb(0.0,1.0,0.0,1.0),type="o")
		if(smooth){
			lines(ParentalScoressmooth,lwd=2,col=rgb(0.0,0.0,0.0,1.0),type="l")
		}
		lines(rep(10.0,length(PopulationScores)),type="l",col=rgb(0.0,1.0,0.0,1.0),lty=2)
	}
	#plot 3
	if(plots[3]){
		plot(InbreadingScores,type="n",...)
		lines(InbreadingScores,lwd=1,col=rgb(0.0,1.0,0.0,0.5))
	}
	if(!plots[4]){
		op <- par(mfrow = c(1,1))
	}
}
