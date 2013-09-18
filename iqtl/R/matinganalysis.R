#####################################################################
#
# matinganalysis.R
#
# copyright parent_check (c) 2009, Danny Arends
# last modified Feb, 2011
# first written Apr, 2009
# 
# Contains: assortedmating
#
######################################################################

assortedmating.init <- function(){
	library(Rexamples)
	workD="D:\\UMCG\\xgap_mating"	
	fileGeno="genotypes_8_matrix.txt"
	fileParents="Parents.csv"
	setwd(workD)
	matrix <- read.table(fileGeno,header=T,row.names=1)
	info <- read.csv(fileParents,header=T)
	assortedmating(matrix,info,1,300,"Geno8")
}

assortedmating.init <- function(){
	workD="D:\\UMCG\\xgap_mating"
	setwd(workD)	
	snps <- as.matrix(read.table("length.txt"))
	ournames <- NULL
	for(x in dir()[grep("marker",dir())]){
		cat("Reading",x,"\n")
		info <- read.table(x,header=T,row.names=1)
		for(searching in 1:nrow(snps)){
			chr <- as.numeric(snps[searching,2])
			if(chr %in% info[,1]){
				bp <- as.numeric(snps[searching,3])
				cat(snps[searching,1:3],"\n")
				smaller <- which(info[,2]<(bp+1000))
				larger <- which(info[,2]>(bp-1000))
				intersect(smaller,larger)
				if(!is.na(intersect(smaller,larger)&&1)){
					cat("->",rownames(info)[intersect(smaller,larger)],"\n")
					ournames <- c(ournames,rownames(info)[intersect(smaller,larger)])
				}
			}
		}
		cat("\n")
	}
	LengthMatrix <- NULL
	for(x in dir()[grep("genotypes",dir())]){
		cat("Reading",x,"\n")
		data <- read.table(x,header=T,row.names=1)
		for(y in ournames){
			if(!is.na(which(rownames(data)==y)&&1)){
			cat("Found ourname:",y,"\n")
				LengthMatrix <- rbind(LengthMatrix,data[which(rownames(data)==y),])
			}
		}
	}
	LengthMatrix
}
	
assortedmating <- function(matrix,info,start=1,number=nrow(matrix)-start,plot=TRUE,verbose=TRUE,logger=FALSE){
	PopulationScores <- NULL		#Teststatistic results holder
	SPopulationScores <- NULL				#Scaled Teststatistic result holder
	InbreadingScores <- NULL				#Inbreading coeficient result holder
	ReturnMatrix <- c("",0,0,0,0)			#Final scoring matrix
	ParentalScores<- NULL					#ParentalMatingsScores
	ParentDirectionS <- NULL				#ParentalMatingSchema LikelyhoodScores
	
	ourcat(logger,"AM","INFO: Analyzing",number,"SNPMarkers\n",verbose=verbose)
	for(locus in start:(start+number)){
		ourcat(logger,"AM","-----------------------------------------------------------------------\n",verbose=verbose)
		ourcat(logger,"AM","Starting with locus:",locus,"\n",verbose=verbose)
		ourcat(logger,"AM","-----------------------------------------------------------------------\n",verbose=verbose)
		p1 <- 0
		p2 <- 0
		p_count <- 0
		c_count <- 0
		Observed <- c(0,0,0)
		Expected <- c(0,0,0)
		Observedtable <- cbind(c(0,0,0),c(0,0,0),c(0,0,0))
		Expectedtable <- cbind(c(0,0,0),c(0,0,0),c(0,0,0))
		teststat1 = 0
		teststat = 0
		ParentDirection <- NULL		#TEMP
		for(ind1 in 1:nrow(info)){
			if(info[ind1,2]==0 && info[ind1,3] !=0){
				#We found a parent
				#cat("Parent found\n")
			}
			if(info[ind1,2] !=0 && info[ind1,3] !=0){
				#Child found
				#cat("Child found\n")
				par1 <- info[ind1,2]
				par2 <- info[ind1,3]
				par1geno <- matrix[locus,par1]
				par2geno <- matrix[locus,par2]
				Observed[(matrix[locus,par1]+1)] <- Observed[(matrix[locus,par1]+1)]+1
				Observed[(matrix[locus,par2]+1)] <- Observed[(matrix[locus,par2]+1)]+1
				p_count<- p_count+2
				Observedtable[(par1geno+1),(par2geno+1)] <- Observedtable[(par1geno+1),(par2geno+1)]+1
				c_count <- c_count+1
			}
			
		}
		n.allells = p_count * 2
		
		ourcat(logger,"AM","Basic info:#Parents:",p_count,",#Children:",c_count,"\n",verbose=verbose)
		ourcat(logger,"AM","Observed: AA(0):\t",Observed[1],"\tAB(1):",Observed[2],"\tBB(2):",Observed[3],"\n",verbose=verbose)
		p1=(2*Observed[1]+Observed[2]) / n.allells
		p2=1-p1
		ourcat(logger,"AM","Calculated: P1:",p1,"\tP2:",p2,"\n",verbose=verbose)
		Expected[1] <- as.integer((p1^2)*p_count)
		Expected[2] <- as.integer((2*p1*p2)*p_count)
		Expected[3] <- as.integer((p2^2)*p_count)
		ourcat(logger,"AM","Expected: AA(0):\t",Expected[1],"\tAB(1):",Expected[2],"\tBB(2):",Expected[3],"\n",verbose=verbose)
		for(x in 1:3){
			if(Expected[x] !=0){
				teststat <- teststat+((Observed[x]-Expected[x])^2/Expected[x])
			}else{
				teststat <- teststat+0
			}
		}
		if(Expected[2]!=0){
			inCoef <- 1-(Observed[2]/Expected[2])
		}else{
			inCoef <- 0
		}
		SPopulationScores <- c(SPopulationScores,teststat*inCoef)
		if(inCoef < 0){
			PopulationScores <- c(PopulationScores,-teststat)
		}else{
			PopulationScores <- c(PopulationScores,teststat)
		}
		InbreadingScores <- c(InbreadingScores,inCoef)
		ourcat(logger,"AM","H/W CHI^2 Test STAT:",teststat,"\n",verbose=verbose)
		if(teststat > 3.84 || teststat < -3.84){
			ourcat(logger,"AM","Population at NOT at equilibrium for locus",locus,"\n",verbose=verbose)
			ourcat(logger,"AM","Inbreading coefficient:",inCoef,"\n",verbose=verbose)
			if(inCoef<0){
				ourcat(logger,"AM"," selection at locus towards HeteroZygote\n",verbose=verbose)
			}else{
				ourcat(logger,"AM"," selection at locus towards HomoZygote\n",verbose=verbose)
			}
		}else{
			ourcat(logger,"AM","Population at equilibrium for locus",locus,"\n",verbose=verbose)
		}
		
		ourcat(logger,"AM","OBSERVED PARENTS\n",verbose=verbose)
		ourcat(logger,"AM","AA\t",Observedtable[1,1],"\t",Observedtable[2,1],"\t",Observedtable[3,1],"\n",verbose=verbose)
		ourcat(logger,"AM","AB\t",Observedtable[1,2],"\t",Observedtable[2,2],"\t",Observedtable[3,2],"\n",verbose=verbose)
		ourcat(logger,"AM","BB\t",Observedtable[1,3],"\t",Observedtable[2,3],"\t",Observedtable[3,3],"\n",verbose=verbose)
		

		ObservedVector <- c(Observedtable[1,1],Observedtable[1,2]+Observedtable[2,1],Observedtable[1,3]+Observedtable[3,1],Observedtable[2,2],Observedtable[3,2]+Observedtable[2,3],Observedtable[3,3])
		for(x in 1:3){
			for(y in 1:3){
				if(x==y){
					Expectedtable[x,y] <- as.integer((Observed[x]/p_count)*(Observed[y]/p_count)*c_count)
				}else{
					Expectedtable[x,y] <- as.integer((Observed[x]/p_count)*(Observed[y]/p_count)*c_count)
				}
			}
		}
		ourcat(logger,"AM","EXPECTED PARENTS\n",verbose=verbose)
		ourcat(logger,"AM","AA\t",Expectedtable[1,1],"\t",Expectedtable[2,1],"\t",Expectedtable[3,1],"\n",verbose=verbose)
		ourcat(logger,"AM","AB\t",Expectedtable[1,2],"\t",Expectedtable[2,2],"\t",Expectedtable[3,2],"\n",verbose=verbose)
		ourcat(logger,"AM","BB\t",Expectedtable[1,3],"\t",Expectedtable[2,3],"\t",Expectedtable[3,3],"\n",verbose=verbose)
		ourcat(logger,"AM","Genotype Parents\t Observed\t Expected\n",verbose=verbose)
		ourcat(logger,"AM","AA*AA\t",Observedtable[1,1],"\t",verbose=verbose)
		ourcat(logger,"AM",Expectedtable[1,1],"\n",verbose=verbose)
		ourcat(logger,"AM","AA*AB\t",Observedtable[1,2]+Observedtable[2,1],"\t",verbose=verbose)
		ourcat(logger,"AM",Expectedtable[1,2]+Expectedtable[2,1],"\n",verbose=verbose)
		ourcat(logger,"AM","AA*BB\t",Observedtable[1,3]+Observedtable[3,1],"\t",verbose=verbose)
		ourcat(logger,"AM",Expectedtable[1,3]+Expectedtable[3,1],"\n",verbose=verbose)
		ourcat(logger,"AM","AB*AB\t",Observedtable[2,2],"\t",verbose=verbose)
		ourcat(logger,"AM",Expectedtable[2,2],"\n",verbose=verbose)
		ourcat(logger,"AM","BB*AB\t",Observedtable[3,2]+Observedtable[2,3],"\t",verbose=verbose)
		ourcat(logger,"AM",Expectedtable[3,2]+Expectedtable[2,3],"\n",verbose=verbose)
		ourcat(logger,"AM","BB*BB\t",Observedtable[3,3],"\t",verbose=verbose)
		ourcat(logger,"AM",Expectedtable[3,3],"\n",verbose=verbose)
		 
		ExpectedVector <- c(Expectedtable[1,1],Expectedtable[1,2]+Expectedtable[2,1],Expectedtable[1,3]+Expectedtable[3,1],Expectedtable[2,2],Expectedtable[3,2]+Expectedtable[2,3],Expectedtable[3,3])
		for(x in 1:6){
			value <- 0
			if(ExpectedVector[x]!=0){
				value <- ((ObservedVector[x]-ExpectedVector[x])^2/ExpectedVector[x])	
			}
			ParentDirection <- c(ParentDirection,value)
			teststat1 <- teststat1+value
		}
		ParentDirection <- c(ParentDirection,(ParentDirection[1]+ParentDirection[4]+ParentDirection[6])/3)
		ParentDirection <- c(ParentDirection,(ParentDirection[2]+ParentDirection[3]+ParentDirection[5])/3)
		ParentDirectionS <- cbind(ParentDirectionS,ParentDirection)
		rownames(ParentDirectionS) <- c("AAxAA","AAxAB","AAxBB","ABxAB","ABxBB","BBxBB","HOMO","HETRO")
		colnames(ParentDirectionS) <- Names("marker",ncol(ParentDirectionS))
		ourcat(logger,"AM","H/W CHI^2 Test STAT:",teststat1,"\n",verbose=verbose)
		if(teststat1 > 7.82){
			ourcat(logger,"AM","ASSORTED MATING suspected (0.05%) on locus",locus,"\n",verbose=verbose)
		}
		ParentalScores <- c(ParentalScores,teststat1)
		ReturnMatrix <- rbind(ReturnMatrix,c(rownames(matrix)[locus],teststat,inCoef,teststat*inCoef,teststat1))
		ourcat(T,"GENO",ParentDirectionS[,(locus-start)],"\n",sep="\t",append=T,verbose=F)
		ourcat(T,"RESULT",ReturnMatrix[(locus-start),],"\n",sep="\t",append=T,verbose=F)
		if(plot){
			temp <- ReturnMatrix
			rownames(temp) <- temp[,1]
			temp <- temp[,-1]
			colnames(temp) <- c("Population_Chi2","Inbreading","Population*Inbreading","Mating_Chi2")
			aa <- list(as.matrix(temp),ParentDirectionS)
			if(nrow(temp) > 10){
				plot.all(aa,smooth=T)
			}else{
				plot.all(aa)
			}
		}
		ourcat(logger,"AM","-----------------------------------------------------------------------\n",verbose=verbose)
		
	}
	rownames(ReturnMatrix) <- ReturnMatrix[,1]
	ReturnMatrix <- ReturnMatrix[,-1]
	ReturnMatrix <- ReturnMatrix[-1,]
	colnames(ReturnMatrix) <- c("Population_Chi2","Inbreading","Population*Inbreading","Mating_Chi2")
	list(ReturnMatrix,ParentDirectionS)
}
