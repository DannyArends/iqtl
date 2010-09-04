\name{plot.genotypes}
\alias{plot.genotypes}
\title{ plot.genotypes - Graphical representation highlighting deviations from expected genotypefreq after randommating }
\description{
	Graphical representation highlighting deviations from expected genotypefreq after randommating
}
\usage{
	plot.genotypes(x, \dots)
}

\arguments{
  \item{x}{ Genotype results of function \link{assortedmating} }
  \item{\dots}{ Additional arguments to plotting function }
}
\details{
	Graphical representation highlighting deviations from expected genotypefreq after randommating
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

}
\seealso{
	\itemize{
		\item \code{\link{assortedmating}} - Assorted mating algorithm
		\item Plot functions:
		\itemize{
			\item \code{\link{plot.popstats}} - Graphical representation of basic population statistics
		}
	}
}
\examples{
	#Create a population at H/W equilibrium
	population <- CreatePopulationHW(300,150)
	
	#Breed assorted (selection between marker:100 and 110) a generation (parents)
	population <- BreedAssorted(population,50,100,110,plot=FALSE)

	#Analyze population using the assortedmating function
	breedingresults <- assortedmating(population$data,population$parentlist[[2]],plot=FALSE,verbose=FALSE)
	
	#To plot genotype selection between parents from breeding results
	plot.genotypes(breedingresults)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
