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
  \item{x}{ Result from the function \link{scoreKinship} }
  \item{start}{ Start at individual }
  \item{num}{ Show this many individuals (DEFAULT: End of the resultset) }
  \item{cutoffsib}{ Sibling cutoff (DEFAULT: value > 3*std-dev) }
  \item{cutoffpar}{ Parental cutoff (DEFAULT: value > 5*std-dev) }
  \item{\dots}{ Additional arguments to plotting function }
}

\details{

}
\value{
	plotting routine, no return value
}
\references{ 

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
		\item \code{\link{scoreKinship}} - Calculates kinshipscores based on SNP markerdata
		\item \code{\link{parents}} - Finds tripplets of individuals with parents from SNP markerdata
		\item \code{\link{kinshipNetwork}} - Created a .SIF network with kinship relations between individuals
	}
}
\examples{
	#Create a population at H/W equilibrium
	population <- CreatePopulationHW(100,300)
	
	#Breed randomly a generation (parents)
	population <- BreedRandom(population,25)
	
	#Score kinship in our breeding example	
	result <- scoreKinship(population$data,plot=FALSE)
	
	#plot the suspected parents
	plot.parentimage(result)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
