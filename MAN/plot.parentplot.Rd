\name{plot.parentplot}
\alias{plot.parentplot}
\title{ plot.parentplot - Graphical representation highlighting possible parents of individuals }
\description{
	Graphical representation highlighting possible parents of individuals 
}
\usage{
	plot.parentplot(x, start = 1, num = (ncol(x)-start),\dots)
}

\arguments{
  \item{x}{ Result from the function \link{scoreKinship} }
  \item{start}{ Start at individual }
  \item{num}{ Show this many individuals }
  \item{\dots}{ Additional arguments to plotting function }
}
\value{
	plotting routine, no return value
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\examples{
	#Create a population at H/W equilibrium
	population <- CreatePopulationHW(100,500)
	
	#Breed randomly a generation (parents)
	population <- BreedRandom(population,50)

	#Score kinship in our breeding example	
	result <- scoreKinship(population$data,plot=FALSE)
	
	plot.parentplot(result,101,4)
	plot.parentplot(result,101,0)
	plot.parentplot(result,103,0)
	plot.parentplot(result,105,0)
}

\keyword{hplot}
