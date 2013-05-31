\name{plotDifCorAtMarker}
\alias{plotDifCorAtMarker}
\title{plotDifCorAtMarker - Plot all the profiles of traits that show differential correlation at a certain marker }
\description{
  Plot all the profiles of traits that show differential correlation at a certain marker
}
\usage{
plotDifCorAtMarker(cross, difCntMatrix, significant=212, lodthreshold=4, marker="YBR008C_211")
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
  \item{difCntMatrix}{ Differential correlation counts matrix. See \code{\link[iqtl]{diffCorAnalysis}} for details.}
  \item{significant}{ Only plot probes which show differential correlation with more then significant other probes }
  \item{lodthreshold}{ Only plot probes which show a lodscore in their QTL profile above lodthreshold }
  \item{marker}{ markername or id of the marker to plot }
}
\details{
  Plot all the profiles of traits that show differential correlation at a certain marker
}
\value{
	QTL mapping results (from \code{\link[qtl]{scanone}}) are returned invisible.
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
