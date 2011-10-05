\name{create.population}
\alias{create.population}
\title{ create.population - Simulation of a population with genetic markers}
\description{
	Simulation of a population genetic markers (not) at H/W equilibrium
}
\usage{
	create.population(individuals,markers,n.geno,normal=TRUE)
}

\arguments{
	\item{individuals}{ Simulate a population with this many individuals}
	\item{markers}{ Simulate this many markers per individuals}
  \item{n.geno}{ Number levels in a single marker genotype}
  \item{normal}{ Should the genotypes follow being at H/W (normal=TRUE) or not (normal=FALSE)}
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\keyword{datagen}
