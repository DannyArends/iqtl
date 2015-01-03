\name{plot.popstats}
\alias{plot.popstats}
\title{plot.popstats - Graphical representation of basic population statistics }
\description{
	Graphical representation of basic population statistics
}
\usage{
	\method{plot}{popstats}(x, smooth=F, num=10, smoothF=median, plots = c(T,T,T,F),\dots)
}

\arguments{
	\item{x}{ Genotype results }
	\item{smooth}{ Draw a smoothes curve ?}
	\item{num}{ Smooth per how many ?}
	\item{smoothF}{ Smoothfunction (eg: mean,median,sum)}
	\item{plots}{ Plot which plots plots <- c(Plot1,Plot2,Plot3,Overplot)}
	\item{\dots}{ Additional arguments to plotting function }
}
\details{
	Plots a graph with basic population statistics
}
\value{
	plotting routine, no return value
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\seealso{
	\itemize{
		\item \code{\link{assortedmating}} - Assorted mating algorithm
		\item Plot functions:
		\itemize{
			\item \code{\link{plot.genotypes}} - Graphical representation of the mating between genotypes per marker
		}
	}
}
\examples{
	#Create a population at H/W equilibrium
	population <- create.population(300,150)
	
	#Breed assorted (selection between marker:100 and 110) a generation (parents)
	#population <- breed.assorted(population,50,100,110,plot=FALSE)

	#Analyze population using the assortedmating function
	#breedingresults <- assortedmating(population$data,population$parentlist[[2]],plot=FALSE,verbose=FALSE)
	
	#To plot basic population statistics
	#plot.popstats(breedingresults,smooth=TRUE,num=10)
}

\keyword{hplot}
