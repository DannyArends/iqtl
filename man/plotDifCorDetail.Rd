\name{plotDifCorDetail}
\alias{plotDifCorDetail}
\title{plotDifCorDetail - Plot a detailed view of a single phenotype in a differential correlation object}
\description{
  plot differential correlation object
}
\usage{
plotDifCorDetail(difCor, pheno.col=1, difCorThreshold=0)
}

\arguments{
  \item{difCor}{ An object of class \code{difCor}. See \code{\link[iqtl]{differentialCorrelation}} for details. }
  \item{pheno.col}{ The phenotype we want to show }
  \item{difCorThreshold}{ Plot only other phenotypes above this threshold }
}
\details{
  Plot a detailed view of a single phenotype in a differential correlation object
}
\value{
	Decreasing ordering of phenotypes relative to the specified pheno.col based on the differential correlation.
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
