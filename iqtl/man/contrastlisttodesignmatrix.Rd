\name{contrastlisttodesignmatrix}
\alias{contrastlisttodesignmatrix}
\title{contrastlisttodesignmatrix - Create a designmatrix from a contrastlist }
\description{
	Create a designmatrix from a contrastlist
}
\usage{
	contrastlisttodesignmatrix(contrastlist, cofactors, verbose=FALSE)
}

\arguments{
  \item{contrastlist}{ List with for each element the marker's contrastmatrix }
  \item{cofactors}{ Cofactors to consider creation of the design matrix }
  \item{verbose}{ Be verbose }
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
