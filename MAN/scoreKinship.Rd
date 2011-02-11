\name{score.kinship}
\alias{score.kinship}
\title{score.kinship - Calculates kinshipscores based on SNP marker data }
\description{
	Calculates kinshipscores based on SNP marker data
}
\usage{
score.kinship(matrix, plot = T, filename = NULL,verbose=FALSE)
}

\arguments{
  \item{matrix}{ matrix[SNPmarker,Individual] }
  \item{plot}{ Should a histogram be plotted (TRUE/FALSE) }
  \item{filename}{ Save result to a file (file=<filename>) }
  \item{verbose}{ Be verbose }
}
\details{
  TODO
}
\value{
	A list with two items:
	\code{Matrix} - Genetic identity
	\code{Matrix} - Genetic distance
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
\seealso{
	\itemize{
		\item \code{\link{find.parents}} - Finds tripplets of individuals with parents from SNP markerdata
		\item \code{\link{kinship.network}} - Created a .SIF network with kinship relations between individuals
		\item Plot functions:
		\itemize{
			\item \code{\link{plot.parentplot}} - Graphical representation highlighting possible parents of individuals
			\item \code{\link{plot.parentimage}} - Graphical representation highlighting possible parents of individuals
		}
	}
}
\examples{
	#Create a population at H/W equilibrium
	population <- create.population(100,300)
	
	#Breed randomly a generation (parents)
	population <- breed.random(population,25)
	
	#Score kinship in our breeding example
	result <- score.kinship(population$data,plot=FALSE)
	
	plot.parentimage(result)
}

\keyword{methods}
