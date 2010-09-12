\name{kinshipNetwork}
\alias{kinshipNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ kinshipNetwork - Created a .SIF network with kinship relations between individuals }
\description{
	Created a .SIF network with kinship relations between individuals 
}
\usage{
kinshipNetwork(result, name="network")
}

\arguments{
  \item{result}{ Result from the function scoreKinship }
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
		\item \code{\link{scoreKinship}} - Calculates kinshipscores based on SNP markerdata
		\item \code{\link{parents}} - Finds tripplets of individuals with parents from SNP markerdata
	}
}
\examples{
	#Create a population at H/W equilibrium
	population <- CreatePopulationHW(300,500)
	
	#Breed randomly a generation (parents)
	population <- BreedRandom(population,50)
	
	#Score kinship in our breeding example
	result <- scoreKinship(population$data,plot=FALSE)
	
	#Print parent 
	pList <- parents(result[[1]])
	pList
	
	#Save a kinship network to a .sif file for use with cytoscape
	kinshipNetwork(result[[1]],"test_network")
}

\keyword{methods}
