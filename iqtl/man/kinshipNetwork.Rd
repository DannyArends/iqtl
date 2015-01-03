\name{kinship.network}
\alias{kinship.network}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ kinship.network - Created a .SIF network with kinship relations between individuals }
\description{
	Created a .SIF network with kinship relations between individuals 
}
\usage{
kinship.network(result, name="network")
}

\arguments{
  \item{result}{ Result from the function score.kinship }
  \item{name}{ Filename for the created networkfile (the extention ".sif" will be added after the name) }
  
}
\details{
	Build a .sif network file containing possible parent <-> child interacytions
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\seealso{
	\itemize{
		\item \code{\link{score.kinship}} - Calculates kinshipscores based on SNP markerdata
		\item \code{\link{find.parents}} - Finds tripplets of individuals with parents from SNP markerdata
	}
}
\examples{
	#Create a population at H/W equilibrium
	population <- create.population(50,100)
	
	#Breed randomly a generation (parents)
	population <- breed.random(population,50)
	
	#Score kinship in our breeding example
	result <- score.kinship(population$data,plot=FALSE)
	
	#Print parent 
	#pList <- parents.network(result[[1]])
	#pList
	
	#Save a kinship network to a .sif file for use with cytoscape
	#kinship.network(result[[1]],"test_network")
}

\keyword{methods}
