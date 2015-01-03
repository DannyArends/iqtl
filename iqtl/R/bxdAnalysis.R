
#setwd("~/Desktop/GeneNetwork")
#bxdGeno       <- read.csv("BXD.geno", skip=6, sep="\t")
#bxdPhenoDescr <- read.csv("BXD.phenome",sep='\t',row.names=1)

getPhenotype <- function(probeset = "13399", DB = "BXDPublish"){
  protocol <- "http"
  url <- "robot.genenetwork.org"
  request <- paste0(protocol, "://", url, "/webqtl/main.py?cmd=get&probeset=", probeset, "&db=", DB, "&format=col")
  return(request, read.csv(, sep = "\t", row.names=1, header=TRUE))
}

possibleBXDstrains <- function(){ return(c("C57BL/6J", "DBA/2J", "B6D2F1", paste0("BXD", 1:103)))}

createPhenotypes <- function(){
  if(!file.exists("BXD.pheno.17sept13")){
    phenotypes <- matrix(NA, length(strains), nrow(bxdPhenoDescr))
    rownames(phenotypes) <- strains
    colnames(phenotypes) <- rownames(bxdPhenoDescr)
    for(x in 1:nrow(bxdPhenoDescr)){
      recordid <- rownames(bxdPhenoDescr)[x]
      phenotype <- getPhenotype(recordid)                                                   # Whole data from geneNetwork
      avail <- which(rownames(phenotype) %in% rownames(phenotypes))                         # We sometimes download data for other BXDs
      phenotypes[rownames(phenotype)[avail], x] <- as.numeric(unlist(phenotype)[avail])     # Storage the available downloaded data
      if(x %% 100 == 0) cat("Downloaded", x,"phenotypes\n")
    }
    write.table(file="BXD.pheno.17sept13", phenotypes,sep='\t',quote=FALSE)
  }else{
    phenotypes <- read.csv("BXD.pheno.17sept13",sep='\t',check.names=FALSE)
  }
  return(phenotypes)
}

# Modify the genotypes
modifyGenotypes <- function(){
	nGeno <- t(bxdGeno[,-c(1:4)])
	nGeno[nGeno=="B"] <- -1
	nGeno[nGeno=="H"] <- 0
	nGeno[nGeno=="D"] <- 1
	nGeno <- apply(nGeno, 2, as.numeric)
	rownames(nGeno) <- colnames(bxdGeno)[-c(1:4)]
	colnames(nGeno) <- bxdGeno[,"Locus"]
}

modifyPhenotypes <- function(){
  noPhenotype <- which(apply(phenotypes,1,function(x){sum(is.na(x))})==ncol(phenotypes))
  phenotypes <- phenotypes[-noPhenotype,]
}

# Map QTLs
mapQTLs <- function(){
  cnt <- 1
  lods <- apply(phenotypes[,1:100], 2, function(phenotype){
	lod <- apply(nGeno, 2, function(marker){
		model <- lm(as.numeric(phenotype[names(marker)]) ~ marker)
		-log10(anova(model)[[5]][1])
	})
	cat(paste0("Done phenotype ", cnt, "\n"));
	cnt <<- cnt+1
	return(lod)
  })
}

toCol <- function(x){
  if(x < 3)  return(rgb(0.8,  0.8, 0.8))
  if(x < 5)  return(rgb(0.7,  0.5, 0.5))
  if(x < 7)  return(rgb(0.7,  0.2, 0.2))
  if(x >= 7) return(rgb(0.7,  0.0, 0.0))
}

getLocsPerChr <- function(bxdGeno, use = c("cM","Mb")){
  use     <- use[1]
  return(lapply(unique(bxdGeno[,1]),function(x){bxdGeno[which(bxdGeno[,1] == x), use]}))
}

getLodsPerChr <- function(bxdGeno, lods, trait = 1){
  return(lapply(unique(bxdGeno[,1]),function(x){lods[which(bxdGeno[,1] == x), trait]}))
}

plotQTL <- function(bxdGeno, lods, trait = 1, use = c("cM","Mb")){
    maplocs <- getLocsPerChr(bxdGeno, use)
    lodlocs <- getLodsPerChr(bxdGeno, lods, trait)
    nchr <- length(unique(bxdGeno[,1]))
    mxx <- max(unlist(lapply(maplocs, max)))
    mnx <- max(unlist(lapply(maplocs, min)))

    plot(c(mnx, mxx),c(0, nchr),t='n', ylab="Chromosome", xlab=paste0("Distance (", use[1],")"))
    chr <- 1
    lapply(maplocs, function(chromosome){
        points(rbind(c(min(chromosome), chr), c(max(chromosome), chr)), t='l', lwd=20, col=rgb(0.9, 0.9, 0.9, 0.9))
        mid <- 1
        lapply(chromosome, function(marker){
            points(rbind(c(marker, (chr-.3)), c(marker, (chr+.3))), t='l',lwd=2, col=toCol(lodlocs[[chr]][mid]))
            mid <<- mid+1
        })
        chr <<- chr +1
    })
    unlist(lapply(lodlocs, function(x){round(max(x), digits = 1)}))
}

createPLots <- function(){
	plotQTL(bxdGeno, lods, 8)

	op <- par(mai=c(3, 4, 0.5, 0.5))
	image(1:ncol(lods),1:ncol(lods), cor(lods), xaxt='n', yaxt='n', ylab="", xlab="", breaks = c(-1, -0.5, -0.2, 0.2, 0.5, 1), col=c("red", "orange", "white", "lightblue", "blue"))
	axis(1, at=1:ncol(lods), colnames(lods), las=2, cex.axis=0.7)
	axis(2, at=1:ncol(lods), colnames(lods), las=1, cex.axis=0.8)
	abline(h=seq(0.5,ncol(lods), 1),col='gray',lty=2)
	abline(v=seq(0.5,ncol(lods), 1),col='gray',lty=2)
}
