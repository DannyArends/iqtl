\name{markercontrasts}
\alias{markercontrasts}
\title{markercontrasts - Calculate per marker contrasts based on the current marker only }
\description{
	Calculate per marker contrasts
}
\usage{
	markercontrasts(genotypes,m=1,verbose=FALSE)
}

\arguments{
	\item{genotypes}{ Genotypes to create contrasts on }
  \item{m}{ Contrast to calculated }
  \item{verbose}{ Verbose output }
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\seealso{
	\itemize{
    \item \code{\link{crosstocontrastlist}} - Create a contrastlist from a cross object
    \item \code{\link{backwardeliminate}} - Create a QTL model by backward elimination of cofactors
    \item \code{\link{contrastqtlmapping}} - QTL mapping routine using multiple linear regression on contrasts
    \item \code{\link{batcheffectcheck}} - Checks for batch effects in the data
		\item \code{\link{batcheffectcorrect}} - Correct batcheffects in data by normalizing group means
	}
}
\examples{
  #TODO
}
\keyword{datagen}
