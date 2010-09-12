\name{scanRF}
\alias{scanRF}
\title{ scanRF - QTL scanning using a single marker anova model }
\description{
	Older version, please see the new mapSM function
}
\usage{
	scanRF(cross=NULL, pheno.col=1,plot=FALSE , compare=FALSE, ...)
}

\arguments{
  \item{cross}{ Cross object to analyse }
  \item{pheno.col}{ one or multiple phenotype columns to analyse }
  \item{plot}{ plot the results }
  \item{compare}{ comapre the results back }
  \item{\dots}{parameters passed to randomfForest }
  }
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  TODO
}
\seealso{
  TODO
}
\examples{
	require(qtl)
	require(randomForest)
	data(listeria)
	data(hyper)	
	aaa <- scanRF(listeria)					#Listeria dataset Phenotype1
	aaa	<- scanRF(hyper)					#Hypertensia dataset Phenotype1
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
