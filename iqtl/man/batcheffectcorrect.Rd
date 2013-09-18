\name{batcheffectcorrect}
\alias{batcheffectcorrect}

\title{ batcheffectcorrect - Correct batch effects in the data }
\description{
	Correct batcheffects in data by normalizing group means using mean centering
}
\usage{
batcheffectcorrect(cross, batchlist, minimumvariance = 0.75)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{batchlist}{ vector with per individual a factor signifying the individuals group }
  \item{minimumvariance}{ minimumvariance threshold, phenotypes with lower variance are excluded }
}
\details{
  TODO
}
\value{
	Matrix with correction values to substract from the phenotype matrix 
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
    \item \code{\link{highvariancephenotypes}} - Select only phenotypes that have a high variance  
		\item \code{\link{batcheffectcheck}} - Checks for batch effects in the data
	}
}
\examples{
	#TODO
}

\keyword{methods}
