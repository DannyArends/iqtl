\name{diffCorPermutation}
\alias{diffCorPermutation}
\title{diffCorPermutation - Calculate differential correlation }
\description{
	Calculate differential correlation by spolitting phenotypes into two groups specified by a genetic marker
}
\usage{
diffCorPermutation(cross, n.perm=10, directory="permutations", verbose=TRUE, \dots)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
  \item{n.perm}{ Number of permutations }
  \item{directory}{ Directory in which to save the output }
  \item{verbose}{ Be verbose }
  \item{\dots}{ Passed to the underlying \code{\link[iqtl]{diffCorAnalysis}} function }
}
\details{
TODO
\itemize{
  \item NOTE: diffCorAnalysis is optimized for multicore machines its uses snow to use two cores
}
}
\value{
difCntPerm object, a list holding difCount matrices from each permutation. To view significance thresholds use 
the \code{\link[iqtl]{significanceThresholds}} function.
}
\references{ 
  TODO
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  TODO
}
\examples{
	#TODO
}

\keyword{methods}
