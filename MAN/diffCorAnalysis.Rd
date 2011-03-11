\name{diffCorAnalysis}
\alias{diffCorAnalysis}
\title{diffCorAnalysis - Calculate differential correlation (and counts) across the genome }
\description{
	Calculate differential correlation (and counts) across the genome
}
\usage{
diffCorAnalysis(cross, minimumvariance=0.01, difCorThreshold=0.4, significant = 5, method="pearson", directory="output", doplot=FALSE, writefile=FALSE, verbose=TRUE)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{minimumvariance}{  Traits with less variance then minimumvariance are removed }
  \item{difCorThreshold}{ Changes in correlation >= difCorThreshold are counted in the difCntMatrix }
  \item{significant}{ Only plot probes which show differential correlation with more then significant other probes }
  \item{method}{ Method of calculating correlation }
  \item{directory}{ Directory in which to save the output }
  \item{doplot}{ Create heatmaps per marker showing differentail correlation across the genotypes }
  \item{writefile}{ Write the differential count matrix file to disk }
  \item{verbose}{ Be verbose }
}
\details{
  TODO
}
\value{
	difCntMatrix object
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
