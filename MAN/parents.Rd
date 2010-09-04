\name{parents}
\alias{parents}

\title{ parents - Finds tripplets of individuals with parents from genetic distance measurements }
\description{
	Finds tripplets of individuals with parents in the output of the \link{scoreKinship} function
}
\usage{
parents(res, file = "test", cutoff = 3)
}

\arguments{
  \item{res}{ Result from the function scoreKinship }
  \item{file}{ Save result to a file (file=<filename>) }
  \item{cutoff}{ cutoff value above which parents are established (DEFAULT: 3*stddev) }
}
\details{

}
\value{
	\code{matrix} - "Individual" "Parent1" "Parent1"
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
		\item \code{\link{scoreKinship}} - Calculates kinshipscores based on SNP markerdata
		\item \code{\link{parents}} - Finds tripplets of individuals with parents from SNP markerdata
		\item \code{\link{kinshipNetwork}} - Created a .SIF network with kinship relations between individuals
	}
}
\examples{
	#Create a population at H/W equilibrium
	population <- CreatePopulationHW(300,500)
	
	#Breed randomly a generation (parents)
	population <- BreedRandom(population,50)
	
	#Score kinship in our breeding example
	result <- scoreKinship(population$data,plot=FALSE)
	
	pList <- parents(result[[1]])
	pList
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
