#
# snow.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: mqmmultitoscanone, mqmmodelsasnetwork
#  Basic scripts for Causal inference
#

inferenetwork <- function(cross,...){
  inferenetwork(traits=pull.pheno(cross),genotypes=pull.geno(cross),...)
}

inferenetwork <- function(traits,genotypes,name="main",QTLresults=NULL,lodqtl=4,loddrop=8){
  traits <- as.matrix(traits)
  genotypes <- as.matrix(genotypes)
  if(is.null(QTLresults)){
    cat("Phase 1: Mapping QTLs\n")
    QTLresults <- mapQTLs(traits,genotypes)
  }else{
    cat("Phase 1: Mapping QTLs... Skipped\n")
  }

  cat("Phase 2: Peak Detection\n")
  peaks <- getpeaks(QTLresults,lodqtl) 
  toinfer <- traitspeaked(peaks)

  cat("Phase 3: Causal Inference\n")
  inferenceprofiles <- calculateconditionallods(traits,genotypes,toinfer)

  cat("Phase 4: Generating network\n")  
  generateSIF(inferenceprofiles,traits,genotypes,type=name,lodqtl,loddrop)
  QTLresults
}


mapQTL <- function(trait, genotypes, markers=seq(1,ncol(genotypes))){
  trait <- as.matrix(trait)
  genotypes <- as.matrix(genotypes)
  LODvalues <- NULL
  for(marker in markers){
    cat(marker,"\n")
    genomodel <- lm(trait ~ as.numeric(genotypes[,marker]))
    Pvalues <- anova(genomodel)[[5]]
    sign <- as.numeric(coefficients(genomodel)[2]/abs(coefficients(genomodel)[2]))
    LODvalues <- c(LODvalues,-log10(Pvalues[1])*sign)
  }
  LODvalues
}

mapQTLs <- function(traits,genotypes){
  traits <- as.matrix(traits)
  qtlprofiles <- NULL
  for(x in 1:ncol(traits)){
    cat(".")
    qtlprofiles <- rbind(qtlprofiles,mapQTL(traits[,x],genotypes))
  }
  cat("\n")
  qtlprofiles
}

inference <- function(traitA, traitB){
  traitA <- as.numeric(traitA)
  traitB <- as.numeric(traitB)
  inferencemodel <- lm(traitA ~ traitB)
  resi <- residuals(inferencemodel)
  resifull <- NULL
  for(x in 1:length(traitA)){
    if(x %in% names(resi)){
      resifull <- c(resifull,resi[which(x == names(resi))])
    }else{
      resifull <- c(resifull,NA)
    }
  }
  resifull
}

getpeaks <- function(qtlprofiles,cutoff = 4.0){
  cat("Starting peak detection above",cutoff,"\n")
  peaks <- vector("list",nrow(qtlprofiles))
  for(x in 1:nrow(qtlprofiles)){
    peak <- FALSE
    curmax <- 0
    curmaxindex <- 1
    marker <- 1
    maximums <- NULL
    for(ab in (qtlprofiles[x,]>cutoff | qtlprofiles[x,]<(-cutoff))){
      if(ab){
        peak <- TRUE
        if(qtlprofiles[x,marker]/abs(qtlprofiles[x,marker]) > 0){
          if(qtlprofiles[x,marker] > curmax){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- marker
          }
        }else{
          if(qtlprofiles[x,marker] < (-curmax)){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- -marker
          }
        }
        if(ncol(qtlprofiles)==marker){
          if(curmax!=0) maximums <- c(maximums,curmaxindex)
        }
      }else{
        if(curmax!=0) maximums <- c(maximums,curmaxindex)
        peak <- FALSE
        curmax <- 0
      }
      marker <- marker+1
    }
    peaks[[x]] <- maximums
  }
  peaks
}

traitspeaked <- function(peaks){
  results <- vector("list",length(unique(unlist(peaks))))
  cnt <- 1
  for(x in unique(unlist(peaks))){
    result <- NULL
    for(y in 1:length(peaks)){
      if(x %in% peaks[[y]]){
        result <- c(result,y)
      }
    }
    results[[cnt]] <- result
    attr(results[[cnt]],"qtl") <- x
    cnt <- cnt + 1
  }
  results
}

calculateconditionallods <- function(traits,genotypes,toinfer){
  inferenceprofiles <- vector("list",length(toinfer))
  for(x in 1:length(toinfer)){
    sign <- attr(toinfer[[x]],"qtl")/abs(attr(toinfer[[x]],"qtl"))
    cat(length(toinfer[[x]]), "traits with a ")
    if(sign < 0 ){
      attr(toinfer[[x]],"qtl") <- (-attr(toinfer[[x]],"qtl"))
      cat("negative")
    }else{
     cat("positive")
    }
    cat(" QTL at marker",attr(toinfer[[x]],"qtl"),"\n")
    inferences <- NULL
    for(traitA in toinfer[[x]]){
      res <- NULL
      for(traitB in toinfer[[x]]){
        if(traitA != traitB){
          res <- c(res,mapQTL(inference(traits[,traitA],traits[,traitB]),genotypes,attr(toinfer[[x]],"qtl")))
        }else{
          res <- c(res,mapQTL(traits[,traitA],genotypes,attr(toinfer[[x]],"qtl")))
        }
      }
      inferences <- rbind(inferences,res)
    }
    rownames(inferences) <- toinfer[[x]]
    colnames(inferences) <- toinfer[[x]]
    inferenceprofiles[[x]] <- inferences
    if(sign < 0){
      attr(toinfer[[x]],"qtl") <- (-attr(toinfer[[x]],"qtl"))
    }
    attr(inferenceprofiles[[x]],"qtl") <- attr(toinfer[[x]],"qtl")
  }
  inferenceprofiles 
}

checkcutoffndrop <- function(prior,past,cutoff,drop){
  if((prior > cutoff | prior<(-cutoff))){
  if(prior > 0){
    if(past < cutoff && abs(past-prior) > drop){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    if(past > (-cutoff) && abs(past-prior) > drop){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  }else{
    return(FALSE)
  }
}

generateSIF <- function(inferenceprofiles,traits,genotypes,type="main",cutoff=2.0,significancedrop=4.0){
  traitnames <- colnames(traits)
  nInd <- 0
  nCau <- 0
  nUnd <- 0
  pEdge <- 0
  nRea <- 0
  cat(file=paste("network",type,".sif",sep=""),"",append=F)
  cat(file=paste("nodes",type,".sif",sep=""),"",append=F)
  for(x in 1:length(inferenceprofiles)){
    traitids <- as.numeric(colnames(inferenceprofiles[[x]]))
    
    for(a in 1:length(traitids)){
      LOD <- inferenceprofiles[[x]][a,a]
      cat(file=paste("nodes",type,".sif",sep=""),colnames(genotypes)[abs(attr(inferenceprofiles[[x]],"qtl"))],"\tMARKER\n",append=TRUE)
      cat(file=paste("nodes",type,".sif",sep=""),traitnames[traitids[a]],"\tTRAIT\n",append=TRUE)
      cat(file=paste("network",type,".sif",sep=""),colnames(genotypes)[abs(attr(inferenceprofiles[[x]],"qtl"))],paste(type,"QTLeffect",sep=""),traitnames[traitids[a]],LOD,"\n",append=TRUE)
      for(b in 1:length(traitids)){
        if(a != b){
          if(checkcutoffndrop(inferenceprofiles[[x]][a,a],inferenceprofiles[[x]][a,b],cutoff,significancedrop)){
            if(checkcutoffndrop(inferenceprofiles[[x]][b,b],inferenceprofiles[[x]][b,a],cutoff,significancedrop)){
              nInd <- nInd+1
            }else{
              cat(file=paste("network",type,".sif",sep=""),traitnames[traitids[b]],paste(type,"Causal",sep=""),traitnames[traitids[a]],abs(inferenceprofiles[[x]][a,a]-inferenceprofiles[[x]][a,b]),"\n",append=TRUE)
              nCau <- nCau+1
            }
          }else{
            if(!checkcutoffndrop(inferenceprofiles[[x]][b,b],inferenceprofiles[[x]][b,a],cutoff,significancedrop)){
              nUnd <- nUnd+1
            }else{
              nRea <- nRea+1
            }
          }
          pEdge <- pEdge+1
        }
      }
    }
  }
  cat("Possible edges: ",pEdge,"\n")
  cat("Causal edges found: ",nCau," (Reactive: ",nRea,")\n")
  cat("Independant edges found: ",nInd,"\n")
  cat("Undecided for ",nUnd," edges\n")
}

