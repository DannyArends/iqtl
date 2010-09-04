\name{scoreKinship}
\alias{scoreKinship}

\title{ scoreKinship - Calculates kinshipscores based on SNPmarker data }
\description{
	Calculates kinshipscores based on SNPmarker data
}
\usage{
scoreKinship(matrix, number = ncol(matrix), plot = T, file = NULL)
}

\arguments{
  \item{matrix}{ matrix[SNPmarker,Individual] }
  \item{number}{ How many individuals of the set should we analyse (DEFAULT: 100, 0 for ALL) }
  \item{plot}{ Should a histogram be plotted (TRUE/FALSE) }
  \item{file}{ Save result to a file (file=<filename>) }
}
\details{

}
\value{
	A list with two items:
	
	\code{Matrix} - Genetic identity
	
	\code{Matrix} - Genetic distance
}
\references{ 

}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{

}
\seealso{
	\itemize{
		\item \code{\link{parents}} - Finds tripplets of individuals with parents from SNP markerdata
		\item \code{\link{kinshipNetwork}} - Created a .SIF network with kinship relations between individuals
		\item Plot functions:
		\itemize{
			\item \code{\link{plot.parentplot}} - Graphical representation highlighting possible parents of individuals
			\item \code{\link{plot.parentimage}} - Graphical representation highlighting possible parents of individuals
		}
	}
}
\examples{
	#Create a population at H/W equilibrium
	population <- CreatePopulationHW(100,300)
	
	#Breed randomly a generation (parents)
	population <- BreedRandom(population,25)
	
	#Score kinship in our breeding example
	result <- scoreKinship(population$data,plot=FALSE)
	
	plot.parentimage(result)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
