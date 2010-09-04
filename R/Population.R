#####################################################################
#
# TestSets.R
#
# copyright parent_check (c) 2009, Danny Arends
# last modified Apr, 2009
# first written Apr, 2009
#
#
######################################################################
CreatePopulationHW <- function(individuals,markers){
	m <- NULL
	pList <- NULL
	m$data <- RnMN(2,individuals,markers)
	m$generation=1
	m$individuals=individuals	
	m$markers=markers
	for(x in 1:individuals){
		pList <- rbind(pList,c(x,0,0))
	}
	m$parentlist <- list(pList)
	colnames(m$data) <- Names("Individual",individuals)
	rownames(m$data) <- Names("Marker",markers)
	m
}

CreatePopulation <- function(individuals,markers){
	m <- NULL
	pList <- NULL
	m$data <- RnM(2,individuals,markers)
	m$generation=1
	m$individuals=individuals	
	m$markers=markers
	for(x in 1:individuals){
		pList <- rbind(pList,c(x,0,0))
	}
	m$parentlist <- list(pList)	
	colnames(m$data) <- Names("Individual",individuals)
	rownames(m$data) <- Names("Marker",markers)
	m
}

make.child <- function(parent1,parent2){
	child <- RnV(2,length(parent1))
	for(x in 1:length(parent1)){
		if(parent1[x]==0 && parent1[x]==parent2[x]){child[x] = parent1[x]}
		if(parent1[x]==2 && parent1[x]==parent2[x]){child[x] = parent1[x]}
		if(parent1[x]==2 && parent2[x]==0){child[x] = 1}
		if(parent1[x]==0 && parent2[x]==2){child[x] = 1}
	}
	child
}

BreedRandom <- function(population,n.child){
	P1 <- RnV(population$individuals-1,n.child)+1
	P2 <- RnV(population$individuals-1,n.child)+1
	result <- NULL
	nextGen <- NULL
	pList <- NULL
	names <- NULL
	for(x in 1:population$individuals){
		pList <- rbind(pList,c(x,0,0))
	}
	for(x in 1:n.child){
		cat("Creating child between:",P1[x],"and",P2[x],"\n")
		pList <- rbind(pList,c(x,P1[x],P2[x]))
		child <- make.child(population$data[,P1[x]],population$data[,P2[x]])
		names <- c(names,paste("Child",population$generation,P1[x],P2[x],sep="."))
		result <- cbind(result,child)
	}
	colnames(result) <- names
	rownames(result) <- Names("Marker",population$markers)
	nextGen$data <- cbind(population$data,result)
	nextGen$individuals <- population$individuals+n.child
	nextGen$markers <- population$marker
	nextGen$generation <- population$generation+1
	nextGen$parentlist <- population$parentlist	
	nextGen$parentlist[[nextGen$generation]] <- pList
	nextGen
}

BreedAssorted <- function(population,n.child,start,stop,...){
	res <- scoreKinship(population$data[start:stop,],population$individuals,...)
	P1 <- RnV(population$individuals-1,n.child)+1
	result <- NULL
	nextGen <- NULL
	pList <- NULL
	names <- NULL

	for(x in 1:n.child){
		P2 <- as.integer(which(max(res[[1]][,P1[x]],na.rm=T)==res[[1]][,P1[x]]))
		cat(x,"Creating child between:",P1[x],"and",P2[1],"\n")
		pList <- rbind(pList,c(x,P1[x],P2[1]))	
		child <- make.child(population$data[,P1[x]],population$data[,P2[1]])
		names <- c(names,paste("Child",population$generation,P1[x],P2[1],sep="."))
		result <- cbind(result,child)		
	}
	colnames(result) <- names
	rownames(result) <- Names("Marker",population$markers)	
	nextGen$data <- cbind(population$data,result)
	nextGen$individuals <- population$individuals+n.child
	nextGen$markers <- population$marker
	nextGen$generation <- population$generation+1
	nextGen$parentlist <- population$parentlist
	nextGen$parentlist[[nextGen$generation]] <- pList
	nextGen
}