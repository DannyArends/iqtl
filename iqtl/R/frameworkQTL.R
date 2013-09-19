#####################################################################
#
# frameworkQTL.R
#
# copyright parent_check (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Mrt, 2009
# 
# Contains: analyseQTL
#
######################################################################


datacheck <- function(data){
	DATAresults <- NULL
	
	DATAresults$individuals <- ncol(data$genotypes)
	DATAresults$markers <- nrow(data$genotypes)
	DATAresults$phenotypes <- length(data$phenotype)

	if(DATAresults$individuals != DATAresults$phenotypes){
		stop("E: individuals != phenotypes")
	}
	
	return(DATAresults)
}

modelQTL <- function(data){
	return(1)
}

mapQTL <- function(data, model){
	mapSM(data,model)
}

permQTL <- function(mapping, data, model, permutations=100){
	permmap <- rep(0,nrow(data$genotypes))
	for(perm in 1:permutations){
		data$phenotype <- sample(data$phenotype)
		map <- mapping(data,model)
		for(marker in 1:nrow(data$genotypes)){
			permmap[marker] <- max(permmap[marker],abs(map[marker]),na.rm=T)
		}
	}
	permmap
}

testdata <- function(ind=10,markers=10){
	data <- NULL
	data$phenotype <- RnV(50,ind)
	data$genotypes <- RnM(2,ind,markers)
	names(data$phenotype) <- Names("IND",ind)
	rownames(data$genotypes) <- Names("M",markers)
	colnames(data$genotypes) <- Names("IND",ind)
	data
}


analyseQTL <- function(data=testdata(),datachecking=datacheck,modelselection=modelQTL,mapping=mapQTL,permutation=permQTL){
	QTLresults <- NULL
	cat("Trait generation\n")
	QTLresults$data     <- data
	cat("Trait check\n")
	QTLresults$check    <- datachecking(QTLresults$data)
	cat("Trait modeling\n")
	QTLresults$model    <- modelselection(QTLresults$data)
	cat("Trait mapping\n")
	QTLresults$traitmap <- mapping(QTLresults$data, QTLresults$model)
	cat("Trait permutation\n")
	QTLresults$permmap  <- permutation(mapping,QTLresults$data, QTLresults$model)
	QTLresults
}
