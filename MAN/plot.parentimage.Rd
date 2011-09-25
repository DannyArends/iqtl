\name{plot.parentimage}
\alias{plot.parentimage}
\title{ plot.parentimage - Graphical representation highlighting possible parents of individuals}
\description{
	Graphical representation highlighting possible parents of individuals
}
\usage{
	plot.parentimage(x, start = 1, num = (ncol(x)-start), cutoffsib = 3, cutoffpar = 5,\dots)
}

\arguments{
  \item{x}{ Result from the function \link{score.kinship} }
  \item{start}{ Start at individual }
  \item{num}{ Show this many individuals (DEFAULT: End of the resultset) }
  \item{cutoffsib}{ Sibling cutoff (DEFAULT: value > 3*std-dev) }
  \item{cutoffpar}{ Parental cutoff (DEFAULT: value > 5*std-dev) }
  \item{\dots}{ Additional arguments to plotting function }
}
\value{
	plotting routine, no return value
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{ 
	num parameter should be larger than 2
}
\seealso{
	\itemize{
		\item \code{\link{score.kinship}} - Calculates kinshipscores based on SNP markerdata
		\item \code{\link{find.parents}} - Finds tripplets of individuals with parents from SNP markerdata
		\item \code{\link{kinship.network}} - Created a .SIF network with kinship relations between individuals
	}
}
\examples{
	#Create a population at H/W equilibrium
	population <- create.population(50, 200)
	
	#Breed randomly a generation (parents)
	population <- breed.random(population,25)
	
	#Score kinship in our breeding example	
	result <- score.kinship(population$data,plot=FALSE)
	
	#plot the suspected parents
	plot.parentimage(result)
}

\keyword{hplot}
