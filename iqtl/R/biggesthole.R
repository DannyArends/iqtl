#####################################################################
#
# biggestholealgorithm.R
#
# copyright (c) 2009, Danny Arends & Joeri v/d Velde
# last modified Sep, 2013
# first written Nov, 2009
# 
# Contains: biggestholealgorithm
#
######################################################################

biggestholealgorithm <- function(x, N = 100){
  if(N < 0) stop("N needs to be larger than 0")
  n <- length(x)
  if(n < 0) stop("x needs to contain items")
  if(n > N) stop("N needs to larger than the number of elements in x")
  bigHole <- 0
  x <- sort(x)
  for(c in 0:n){
    if(c==0){
	  g1 <- 0
	}else{
	  g1 <- x[c]
	}
	if(c==n){
	  g2 <- N
	}else{
	  g2 <- x[c+1]
	}
	if(g1 < 0 || g1 > N) stop("Element",c,"value",g1,"outside domain")
	if(g2 < 0 || g2 > N) stop("Element",c,"value",g2,"outside domain")
	if(bigHole < (g2-g1)){ bigHole <- (g2-g1) }
  }
  cat("Maximum expected hole:", round(N/n), "\n")
  cat("Observed largest hole:", bigHole, "\n")
  return (round(N/n)/(bigHole))
}


