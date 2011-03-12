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
\itemize{
  \item NOTE: diffCorAnalysis is optimized for multicore machines its uses snow to use two cores
}
}
\value{
difCor object, a list holding:
\itemize{
  \item [[1]] Differential correlation matrix (trait x trait) based on the differences in correlation^2
  \item [[2]] Correlation matrix (trait x trait) based only on individuals carrying genotype A at the specified marker
  \item [[3]] Correlation matrix (trait x trait) based only on individuals carrying genotype B at the specified marker
  \item [[4]] A vector of counted scores for each phenotype above the difCorThreshold
}
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
