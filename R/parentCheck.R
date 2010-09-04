#####################################################################
#
# parentCheck.R
#
# copyright parent_check (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Mrt, 2009
# 
# Contains: parentCheck
# Input: A matrix with in the columns individuals and in the rows markers with genotypes for each individual
# Returns: A matrix with scores indicating kinship between individuals
# Output: a netwerk.sif file with the generated kinship network
#
# Contains: parentCheck2
# Input: A matrix with in the columns individuals and in the rows markers with genotypes for each individual
# Returns: A matrix with scores indicating kinship between individuals
# Output: a netwerk.sif file with the generated kinship network
# 
# Contains: kinshipNetwork
# Input: A matrix with scores indicating kinship between individuals
# Returns: Nothing
# Output: A netwerk_new.sif file with the generated kinship network
#
# Furthermore: find.parent,parents
#
######################################################################

scoreKinship <- function(matrix,number=ncol(matrix),plot=T,file=NULL){
	result1 <- NULL
	result2 <- NULL
	b <- NULL
	
	for(x in 1:number){
		result1 <- cbind(result1,rep(NA,number))
		result2 <- cbind(result2,rep(NA,number))
	}
	colnames(result1) <- colnames(matrix)[1:number]
	rownames(result1) <- colnames(matrix)[1:number]
	colnames(result2) <- colnames(matrix)[1:number]
	rownames(result2) <- colnames(matrix)[1:number]
	
	pb <- txtProgressBar(min=1,max=number,style=3)
	for(ind1 in 1:number){
		setTxtProgressBar(pb, ind1)
		start <- proc.time()
		others <- 1:number
		others <- others[-ind1]
		for(ind2 in others){
			score1 <- sum(matrix[,ind1]==matrix[,ind2])
			score1 <- score1/nrow(matrix)
			score2 <- sum(abs(matrix[,ind1]-matrix[,ind2]))
			score2 <- score2/(nrow(matrix)*2)			
			result1[ind1,ind2] <- score1
			result2[ind1,ind2] <- score2
		}
		if(plot){
			plot.histkinship(list(result1,result2),ind1)
		}
	}
	cat("\n")	#FIX: ProgressBar doesn't print an endline	
	if(!is.null(file)){
		cat("INFO: Writing score files\n")
		write.csv(result1,file=paste("R1",file,sep=""))
		write.csv(result2,file=paste("R2",file,sep=""))
	}
	list(result1,result2)
}

parentCheck <- function(matrix,number=100){
	result <- NULL
	cat(paste(""), file = "netwerk.sif",append=FALSE)
	for(x in 1:number){result <- cbind(result,rep(NA,number))}
	colnames(result) <- paste("Ouder",1:ncol(result),sep="")
	rownames(result) <- paste("Kind",1:nrow(result),sep="")
	pb <- txtProgressBar(min=1,max=ncol(result),style=3)
	for(kind in 1:number){
		setTxtProgressBar(pb, kind)
		start <- proc.time()
		ouders <- 1:number
		ouders <- ouders[-kind]
		res <- NULL
		m_score_1 <- -1
		m_ouder_1 <- -1
		m_score_2 <- -1
		m_ouder_2 <- -1
		for(ouder in ouders){
			score <- sum((matrix[,kind]==2)==(matrix[,ouder]==2)) + sum((matrix[,kind]==0)==(matrix[,ouder]==0))
			score = score - (sum((matrix[,kind]==0)==(matrix[,ouder]==2)) + sum((matrix[,kind]==2)==(matrix[,ouder]==0)))
			result[kind,ouder] <- score
			if(score > m_score_1){
				m_score_1 <- score
				m_ouder_1 <- ouder
			}else{
				if(score > m_score_2){
					m_score_2 <- score
					m_ouder_2 <- ouder
				}
			}
		}
		hist(result[kind,],main=paste("Scores voor ouders, kind=",kind),xlab="Parental scores")
		meani <- mean(result[kind,],na.rm=T)
		sdi <- sd(result[kind,],na.rm=T)
		cutoff <- meani+3*sdi
		#cat("--------------------------------------------------------------------------\n")
		#cat("Individu:",kind,"\n")
		#if(m_score_1 > cutoff && m_score_2 > cutoff){
		#	cat("Is waarschijnlijk een kind van ouders:\n")
		#	cat(paste(kind,"->",m_ouder_1,"\n"), file = "netwerk.sif",append=TRUE)
		#	cat(paste(kind,"->",m_ouder_2,"\n"), file = "netwerk.sif",append=TRUE)
		#	cat("Hoogste verwantschap 1:",m_ouder_1,"met:",m_score_1,"\n")
		#	cat("Hoogste verwantschap 2:",m_ouder_2,"met:",m_score_2,"\n")
		#}else{
		#	if(m_score_1 < cutoff && m_score_2 < cutoff){
		#		cat("Waarschijnlijk geen verwantschap in set\n")
		#	}
		#	if(m_score_1 > cutoff){
		#		cat("Is waarschijnlijk een ouder:\n")
		#		cat(paste(kind,"->",m_ouder_1,"\n"), file = "netwerk.sif",append=TRUE)
		#		cat("Hoogste verwantschap:",m_ouder_1,"met:",m_score_1,"\n")
		#	}
		#	if(m_score_2 > cutoff){
		#		cat("Is waarschijnlijk een ouder:\n")
		#		cat(paste(kind,"->",m_ouder_2,"\n"), file = "netwerk.sif",append=TRUE)
		#		cat("Hoogste verwantschap:",m_ouder_2,"met:",m_score_2,"\n")
		#	}	
		#}
		#end <- proc.time()
		#cat("--------------------------------------------------------------------------\n")
		#cat("INFO: Calculation of x=",kind," ",round((end-start)[3], digits=3)," seconds\n",sep="")
		#cat("--------------------------------------------------------------------------\n")
	}
	cat("\n")
	result
}

kinshipNetwork <- function(result,name="network"){
	cat(paste(""), file = paste(name,".sif",sep=""),append=FALSE)
	pb <- txtProgressBar(min=1,max=ncol(result),style=3)
	number <- ncol(result)
	for(kind in 1:number){
		setTxtProgressBar(pb, kind)
		start <- proc.time()
		done <- NULL
		meani <- mean(result[kind,],na.rm=T)
		sdi <- sd(result[kind,],na.rm=T)
		for(x in 10:4){
			cutoff_v <- meani+x*sdi
			a <- which(result[kind,]>=cutoff_v)
			for(y in a){
				if(y %in% done){
					#cat("Interaction already included\n")
				}else{
					cat(paste(kind,paste("i",x,sep=""),y,"\n"), file = paste(name,".sif",sep=""),append=TRUE)
					done <- c(done,y)
				}
			}
		}
	}
	cat("\n")	
}

find.parent <- function(res,individual,cutoff){
	meani <- mean(res[individual,],na.rm=T)
	sdi <- mean(sd(res[individual,],na.rm=T))
	cutoff_p = meani+cutoff*sdi
	r <- NULL
	a <- which(res[individual,] >= cutoff_p)
	#cat("IND:",individual,"Num above threshold (",cutoff_p,"):",length(a),"\n")
	r <- c(individual,0,0)
	if(length(a)==1){
		#cat("IND:",individual,"Num above threshold (",cutoff_p,"):",length(a),"\n")
		r <- c(individual,0,a[1])
	}
	if(length(a)==2){
		r <- c(individual,a[1],a[2])
	}
	r
}

parents <- function(res,file="test",cutoff=3){
	r <- NULL
	res <- as.matrix(res)
	pb <- txtProgressBar(min=1,max=ncol(res),style=3)
	for(x in 1:dim(res)[1]){
		setTxtProgressBar(pb, x)
		r <- rbind(r,find.parent(res,x,cutoff))
	}
	rownames(r) <- r[,1]
	r <- r[,-1]
	colnames(r) <- c("Ouder1","Ouder2")
	write.csv(r,file=paste(file,".csv",sep=""))
	cat("\n")
	cat("INFO: Number of individuals that got parents assigned to them:",(length(r)/2)-sum(r[,1]==0),"\n")
	for(x in 1:dim(r)[1]){
		if(r[x,1]!=0){
			o1 <- r[x,1]
			o2 <- r[x,2]
			#cat("INFO: IND:",x,"With parents:",o1,o2,"\n")
			if(r[o1,2]==x && r[o1,1]==0 && r[o2,2]==x && r[o2,1]==0){
				#cat("INFO: HOLDS\n")	
				cat(rownames(r)[x],"->",r[x,1],"\n",file=paste(file,".sif",sep=""),append=T)
				cat(rownames(r)[x],"->",r[x,2],"\n",file=paste(file,".sif",sep=""),append=T)
			}else{
				#cat("INFO: Failed second test")
			}
		}
	}
	r
}



doDIR <- function(markers=1000){
	pb <- txtProgressBar(min=1,max=markers,style=3)
	for(x in grep("genotypes",dir())){
		cat("INFO: Opening File:",dir()[x],"\n")
		data <- read.table(dir()[x],header=T,row.names=1)
		chosen <- NULL
		while(length(chosen) < markers){
			add <- as.integer(runif(1, min=1, max=nrow(data)))
			if(!add %in% chosen){
				setTxtProgressBar(pb, length(chosen))
				chosen <-c(add,chosen)
			}
		}
		cat("\n")	#FIX: ProgressBar doesn't print an endline	
		cat("INFO: Gonna score kinship\n")
		scoreKinship(data[chosen,],ncol(data),T,file=dir()[x])	
		cat("INFO: File:",dir()[x],"done\n")		
	}
}

avgResults <- function(){
	setwd("D:\\GBIC\\Yoeri&Morris\\parentanalyse")
	result <- NULL
	data <- NULL
	for(x in grep("R2genotypes",dir())){	
		cat("INFO: Opening File:",dir()[x],"\n")
		data <- read.csv(dir()[x],header=T,row.names=1)
		if(x==25){
			result <- data
		}else{
			result <- result+data
		}
	}
	avg_res <- result/length(grep("R2genotypes",dir()))
	write.csv(avg_res,)
}

loadAvgResults <- function(){
	setwd("D:\\GBIC\\Yoeri&Morris\\parentanalyse")
	data <- read.csv("R1avg.csv",header=T,row.names=1)
	data
}
