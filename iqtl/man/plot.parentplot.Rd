\name{plot.parentplot}
\alias{plot.parentplot}
\title{plot.parentplot - Graphical representation highlighting possible parents of individuals }
\description{
	Graphical representation highlighting possible parents of individuals 
}
\usage{
	\method{plot}{parentplot}(x, start = 1, num = (ncol(x)-start),\dots)
}

\arguments{
  \item{x}{ Result from the function \link{score.kinship} }
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
	population <- create.population(100,300)
	
	#Breed randomly a generation (parents)
	population <- breed.random(population,50)

	#Score kinship in our breeding example	
	result <- score.kinship(population$data,plot=FALSE)
	
	plot.parentplot(result,101,4)
	plot.parentplot(result,101,0)
	plot.parentplot(result,103,0)
	plot.parentplot(result,105,0)
}

\keyword{hplot}
