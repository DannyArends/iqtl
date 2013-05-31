\name{eQCLpermute}
\alias{eQCLpermute}
\title{eQCLpermute - permutation routine for eQCL analysis }
\description{
	permutation routine for eQCL analysis
}
\usage{
eQCLpermute(cross, n.perm=10, directory="permutations", verbose=TRUE, \dots)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
  \item{n.perm}{ Number of permutations to perform (recommended >1000) }
  \item{directory}{ Directory in which to save the output }
  \item{verbose}{ Be verbose }
  \item{\dots}{ Passed to the underlying \code{\link[iqtl]{eQCLscan}} function }
}
\details{
TODO
\itemize{
  \item NOTE: The underlying function eQCLscan is optimized for multicore machines its uses the snow package
}
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
  require(iqtl)
  data(multitrait)
  multi_perm <- eQTLpermute(multitrait)
  significance.eQCL(multi_perm)
}

\keyword{methods}
