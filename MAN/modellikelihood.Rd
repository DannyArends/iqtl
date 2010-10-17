\name{modellikelihood}
\alias{modellikelihood}
\title{modellikelihood - Get the unadjusted log likelyhood of a certain model }
\description{
	Get the unadjusted log likelyhood of a certain model using multiple linear regression
}
\usage{
	modellikelihood(designmatrix,y,weight=rep(1,nrow(designmatrix)),verbose=FALSE)
}

\arguments{
  \item{designmatrix}{ Designmatrix with samples in rows, and variables in columns }
  \item{y}{ Response vector (length: nrow) }
  \item{weight}{ Weight of the samples (length: nrow) }
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
    \item \code{\link{contrastqtlsignificance}} - QTL permutation routine using multiple linear regression on contrasts
	}
}
\examples{
  #TODO
}
\keyword{methods}
