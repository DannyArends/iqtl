\name{plotDifCntProfile}
\alias{plotDifCntProfile}
\title{plotDifCntProfile - plot differential correlation count profile for a single probe}
\description{
  plot differential correlation count profile for a given probe
}
\usage{
plotDifCntProfile(cross, difCntMatrix, pheno.col=1, significanceThresholds=NULL, addQTL=FALSE, \dots}
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{difCntMatrix}{ Differential correlation counts matrix. See \code{\link[iqtl]{diffCorAnalysis}} for details.}
  \item{pheno.col}{ Which phenotype do we want to show }
  \item{significanceThresholds}{ The resulting object from significanceThresholds to a DifPermObject }
  \item{addQTL}{ Add a QTL profile to the plot (scanone is used to scan for QTL) }
  \item{\dots}{ Passed to the underlying plot function }
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
