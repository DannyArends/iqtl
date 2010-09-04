# pval_lod_converter.r
#
# copyright (c) Danny Arends, Frank Johannes, Joeri van der Velde 2009
#
# last modified may,2009
# first written may, 2009
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
 
pvaltolod <- function(pvalue,accuracy=(1/100),df=1,lower.tail=FALSE){
	t <- seq(0,100,accuracy)
	a <- pchisq(t,df=df,lower.tail=lower.tail)
	res <- NULL
	for(i in 1:length(a)){
		if(!(a[i] > pvalue) && is.null(res)){
		 res <- (i*accuracy)-accuracy
		 res <- log(exp(res/2),10)
		}
	}
	res
}

lodtopval <- function(lod,df=1,lower.tail=FALSE){
	likelihoodratio<-10^abs(lod)
	chisquarevalue<-2*log(likelihoodratio)
	pvalue <- pchisq(chisquarevalue, df=df, lower.tail=lower.tail)
	pvalue
}

examples <- function(){
	# Pvalue cutoffs
	eqtl = 5.547315e-05
	nmr = 0.0005353726
	lcms= 0.0003983748

	# convert Pvalue to LOD score
	eqtl_lod <- pvaltolod(eqtl)
	eqtl_lod_np <- pvaltolod(eqtl,accuracy=(1/5))
	eqtl_lod_nnp <- pvaltolod(eqtl,accuracy=(10/1))	
	nmr_lod <- pvaltolod(nmr)
	lcms_lod <- pvaltolod(lcms)

	# and back: convert LOD score to Pvalue
	eqtl_pval <- lodtopval(eqtl_lod)
	eqtl_pval_np <- lodtopval(eqtl_lod_np)
	eqtl_pval_nnp <- lodtopval(eqtl_lod_nnp)	
	nmr_pval <- lodtopval(nmr_lod)
	lcms_pval <- lodtopval(lcms_lod)
	cat("Name\tAccuracy\tPval (p)\tLOD (c)\t\tPval (c)\n")
	cat("eqtl\t1/100\t",eqtl,eqtl_lod,eqtl_pval,"\n",sep="\t")
	cat("eqtl\t1/5\t",eqtl,eqtl_lod_np,eqtl_pval_np,"\n",sep="\t")
	cat("eqtl\t10/1\t",eqtl,eqtl_lod_nnp,eqtl_pval_nnp,"\n",sep="\t")
}
