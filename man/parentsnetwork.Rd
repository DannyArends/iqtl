\name{parents.network}
\alias{parents.network}
\title{parents.network - Finds tripplets of individuals with parents from genetic distance measurements }
\description{
	Finds tripplets of individuals with parents in the output of the \link{score.kinship} function
}
\usage{
parents.network(res, file = "test", cutoff = 3)
}

\arguments{
  \item{res}{ Result from the function score.kinship }
  \item{file}{ Save result to a file (file=<filename>) }
  \item{cutoff}{ cutoff value above which parents are established (DEFAULT: 3*stddev) }
}
\value{
	\code{matrix} - "Individual" "Parent1" "Parent1"
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
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
	population <- create.population(50,100)
	
	#Breed randomly a generation (parents)
	population <- breed.random(population,50)
	
	#Score kinship in our breeding example
	result <- score.kinship(population$data,plot=FALSE)
	
	pList <- parents.network(result[[1]])
	pList
}

\keyword{methods}
