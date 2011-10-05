\name{autonormalize}
\alias{autonormalize}
\title{autonormalize - Automated per phenotype normalization }
\description{
	Automated per phenotype normalization strategy trying probit, logit, sqrt, reciproce
}
\usage{
	autonormalize(cross,transformations=NULL,verbose=TRUE)
}

\arguments{
	\item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
	\item{transformations}{ Vector with normalization methods for each phenotype }
	\item{verbose}{ Verbose output }
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\seealso{
	\itemize{
    \item \code{\link{batcheffectcheck}} - Checks for batch effects in the data
		\item \code{\link{batcheffectcorrect}} - Correct batcheffects in data by normalizing group means
    \item \code{\link{backwardeliminate}} - Create a QTL model by backward elimination of cofactors
    \item \code{\link{contrastqtlmapping}} - QTL mapping routine using multiple linear regression on contrasts
	}
}
\examples{
  data(multitrait)
  multitrait$pheno <- data.frame(matrix(
  c(rnorm(810),
    runif(810),
    exp(rnorm(810)),
    sqrt(1+runif(810))
    ),162,20))
  results <- autonormalize(multitrait)
}
\keyword{datagen}
