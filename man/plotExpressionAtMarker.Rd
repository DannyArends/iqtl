\name{plotExpressionAtMarker}
\alias{plotExpressionAtMarker}
\title{plotExpressionAtMarker - Box plot of expressions split by genotype }
\description{
  Plots the expression specified by pheno.col but it splits it based on the genotype at marker
}
\usage{
plotExpressionAtMarker(cross, pheno.col=1, marker="YBR008C_211")
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
  \item{pheno.col}{ Phenotype column to plot }
  \item{marker}{ Marker used to split the individuals based on genotype }
}
\details{
  Plots the expression specified by pheno.col but it splits it based on the genotype at marker
}
\value{
	List with expressions split in both genotype groups
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
