\name{scaledownPhenotypes}
\alias{scaledownPhenotypes}
\title{scaledownPhenotypes - Remove phenotypes with low variation }
\description{
	Calculates kinshipscores based on SNP marker data
}
\usage{
scaledownPhenotypes(cross, minimumvariance = 0.1,verbose=FALSE)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{minimumvariance}{ Traits with less variance then minimumvariance are removed }
  \item{verbose}{ Be verbose }
}
\details{
  TODO
}
\value{
	An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details containing only phenotypes with variance >= minimumvariance
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
