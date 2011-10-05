\name{plotDifCor}
\alias{plotDifCor}
\title{plotDifCor - plot differential correlation object}
\description{
  plot differential correlation object
}
\usage{
plotDifCor(difCor, difCorThreshold=0.5, significant = 0, \dots)
}

\arguments{
  \item{difCor}{ An object of class \code{difCor}. See
                \code{\link[iqtl]{differentialCorrelation}} for details. }
  \item{difCorThreshold}{ Changes in correlation >= difCorThreshold are counted in the difCntMatrix }
  \item{significant}{ Only plot probes which show differential correlation with more then significant other probes }
  \item{\dots}{ Arguments passed to image }
}
\details{
  TODO
}
\value{
	difCntMatrix
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
