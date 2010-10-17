\name{autonormalize}
\alias{autonormalize}
\title{ autonormalize - Automated per phenotype normalization }
\description{
	Automated per phenotype normalization strategy trying probit, logit, sqrt, reciproce
}
\usage{
	autonormalize(cross,transformations=NULL,verbose=TRUE){
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
\examples{
  #TODO
}
\keyword{datagen}
