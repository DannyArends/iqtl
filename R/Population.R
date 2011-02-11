#####################################################################
#
# population.R
#
# copyright population (c) 2009, Danny Arends
# last modified Feb, 2011
# first written Apr, 2009
#
# 
# Contains: random.genotypes, create.population, create.child, breed.random, breed.assorted
#
######################################################################

random.genotypes <- function(individuals, markers, n.geno=2, normal=TRUE){
	if(normal){
		aa <- RnMN(n.geno,individuals,markers)
	}else{
		aa <- RnM(n.geno,individuals,markers)
	}
	rownames(aa) <- Names("marker",markers)
	colnames(aa) <- Names("ind",individuals)
	aa
}

create.population <- function(individuals, markers, n.geno=2, normal=TRUE){
	m <- NULL
	parentlist <- NULL
	m$data <- random.genotypes(individuals,markers,n.geno,normal)
	m$generation=1
	m$individuals=individuals	
	m$markers=markers
	for(x in 1:individuals){
		parentlist <- rbind(parentlist,c(x,0,0))
	}
	m$parentlist <- list(parentlist)
	colnames(m$data) <- Names("Individual",individuals)
	rownames(m$data) <- Names("Marker",markers)
	m
}

create.child <- function(parent1,parent2){
	child <- RnV(2,length(parent1))
	for(x in 1:length(parent1)){
		if(parent1[x]==0 && parent1[x]==parent2[x]){child[x] = parent1[x]}
		if(parent1[x]==2 && parent1[x]==parent2[x]){child[x] = parent1[x]}
		if(parent1[x]==2 && parent2[x]==0){child[x] = 1}
		if(parent1[x]==0 && parent2[x]==2){child[x] = 1}
	}
	child
}

breed.random <- function(population,n.child){
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
		child <- create.child(population$data[,P1[x]],population$data[,P2[x]])
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

breed.assorted <- function(population, n.child, start, stop, ...){
	res <- score.kinship(population$data[start:stop,],population$individuals,...)
	P1 <- RnV(population$individuals-1,n.child)+1
	result <- NULL
	nextGen <- NULL
	pList <- NULL
	names <- NULL

	for(x in 1:n.child){
		P2 <- as.integer(which(max(res[[1]][,P1[x]],na.rm=T)==res[[1]][,P1[x]]))
		cat(x,"Creating child between:",P1[x],"and",P2[1],"\n")
		pList <- rbind(pList,c(x,P1[x],P2[1]))	
		child <- create.child(population$data[,P1[x]],population$data[,P2[1]])
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
