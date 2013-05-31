\name{load.eQCLpermute}
\alias{load.eQCLpermute}
\title{load.eQCLpermute - permutation routine for eQCL analysis }
\description{
  Reload an eQCL permutation analysis from disk
}
\usage{
load.eQCLpermute(cross, n.perm=10, directory="permutations", verbose=TRUE, \dots)
}

\arguments{
  \item{directory}{ Directory in which to save the output }
}
\details{
TODO
}
\value{
eQTLperm object, a list holding eQTL matrices from each permutation. To view significance thresholds use 
the \code{\link[iqtl]{significance.eQCL}} function.
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
  #multi_perm <- load.eQTLpermute("directory")
  #significance.eQCL(multi_perm)
}

\keyword{methods}
