\name{differentialCorrelation}
\alias{differentialCorrelation}
\title{differentialCorrelation - Calculate differential correlation }
\description{
	Calculate differential correlation by spolitting phenotypes into two groups specified by a genetic marker
}
\usage{
differentialCorrelation(cross, marker, difCorThreshold=0.5, method="pearson", directory="output", saveRdata=FALSE)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{marker}{ Which marker should we analyse }
  \item{difCorThreshold}{ Changes in correlation >= difCorThreshold are counted in the difCntMatrix }
  \item{method}{ Method of calculating correlation }
  \item{directory}{ Directory in which to save the output }
  \item{saveRdata}{ Save the correlation matrices as Rdata object to the directory }
}
\details{
  TODO
}
\value{
	difCor object, a list holding:
  [[1]] Differential correlation based on genotype at this marker
  [[2]] correlation group 1 based on genotype at this marker
  [[3]] correlation group 2 based on genotype at this marker
  [[4]] difCountVector based on genotype at this marker
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
