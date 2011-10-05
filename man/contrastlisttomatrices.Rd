\name{contrastlisttomatrices}
\alias{contrastlisttomatrices}
\title{contrastlisttomatrices - Create a designmatrix and nullmatrix from a contrastlist }
\description{
	Create a designmatrix and nullmatrix from a contrastlist
}
\usage{
	contrastlisttomatrices(contrastlist, m, cofactors = NULL, verbose=FALSE)
}

\arguments{
	\item{contrastlist}{ List with for each element the marker's contrastmatrix }
  \item{m}{ Contrast that is is being nulled for the nullmodel }
  \item{cofactors}{ Cofactors to consider creation of the design matrix }
  \item{verbose}{ Verbose output }
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\seealso{
	\itemize{
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
