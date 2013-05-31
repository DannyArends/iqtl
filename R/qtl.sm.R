#####################################################################
#
# smQTL.R
#
# copyright (c) 2009, Danny Arends
# last modified March, 2008
# first written Unknown
#
# Contains: scanSM
#
######################################################################
#   Description: Single marker analysis for QTL mapping
#   Arguments:
#       DATA
#		    data$genotypes:	Matrix of genotype. Rows represent markers and the columns represent individuals
#       data$traits: 	Matrix of phenotype. Rows represent traits and the columns represent invididuals
#       file: 	The output filename
#   Details:
#       Does an anova test to calculate LOD scores for linkage likelihood of phenotype values and marker genotypes at each marker position
#   Output:
#       Matrix of -logp values for linkage likelihood. the rows represent markers colomns represent traits.
#       Positive/negative sign is added to represent the direction of additive effect
#       positive: phenotype values 2 > phenotype values genotype 1.
######################################################################

mapSM <- function(data,model){
	genotypes <- data$genotypes
	traits <- data$phenotype
	markers <- nrow(genotypes)
	map <- NULL
	for (marker in 1:markers){
		model <- lm(traits~genotypes[marker,])
		myanova <- anova(model)[[5]][1]
		if (model[[1]][2]<0){
			ef <- -1
		}else{
			ef <- 1
		}
		map <- c(map, -log10(myanova)*ef)
	}
	names(map) <- Names("M",markers)
	map
}


scanSM <- function(genotypes,traits,file=NULL, sep=","){
	n.marker <- nrow(genotypes)
	if (is.vector(traits)){
		traits <- t(as.matrix(traits))
		n.traits <- nrow(traits)
	}
	if (!is.null(rownames(traits))){
		name.traits <- rownames(traits)
	}else{
		name.traits <- paste("trait", 1:n.traits, sep="")
	}
	if (!is.null(rownames(genotypes))){
		name.marker <- rownames(genotypes)
	}else{
		name.marker <- paste("M", 1:n.marker, sep="")
	}
	if(!is.null(file)){
		cat(rownames(genotypes), "\n",file=file, sep=sep, append=F)
	}
	lod <- NULL
	for (i in 1:n.traits){
		qtl.p<-NULL
		for (j in 1:n.marker){
			model <- lm(traits[i,]~genotypes[j,])
			myanova <- anova(model)[[5]][1]
			if (model[[1]][2]<0){
				ef <- -1
			}else{
				ef <- 1
			}
			qtl.p <- c(qtl.p, -log10(myanova)*ef)
		}
		if(!is.null(file)){
			cat (name.traits[i], qtl.p,"\n", file=file, sep=sep, append=T)
		}
		lod <- rbind(lod, qtl.p)
	}
	dimnames(lod)<-list(name.traits, name.marker)
	lod
}
