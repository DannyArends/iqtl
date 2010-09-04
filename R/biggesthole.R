#####################################################################
#
# biggestholealgorithm.R
#
# copyright (c) 2009, Danny Arends & Joeri v/d Velde
# last modified Nov, 2009
# first written Nov, 2009
# 
# Contains: biggestholealgorithm
#
######################################################################

biggestholealgorithm <- function(N = 100, list){
	if(N < 0) stop("N needs to be larger than 0")
	n <- length(list)
	if(n < 0) stop("list needs to contain items")
	if(n > N) stop("N needs to larger than the number of elements in list")
	hetgrootstegat <- 0
	list <- sort(list)
	for(x in 0:n){
		if(x==0){
			g1 <- 0
		}else{
			g1 <- list[x]
		}
		if(x==n){
			g2 <- N
		}else{
			g2 <- list[x+1]
		}
		if(g1 < 0 || g1 > N) stop("list element",x,"value",g1,"outside domain")
		if(g2 < 0 || g2 > N) stop("list element",x,"value",g2,"outside domain")
		if(hetgrootstegat < (g2-g1)){
			hetgrootstegat <- (g2-g1)
		}
	}
	cat("Maximum expected hole:",round(N/n),"\n")
	cat("Observed largest hole:",hetgrootstegat,"\n")
	if(hetgrootstegat)
	return (round(N/n)/(hetgrootstegat))
}


