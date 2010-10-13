\name{batcheffectcheck}
\alias{batcheffectcheck}

\title{ batcheffectcheck - Checks for batch effects in the data }
\description{
	Checks for batch effects in the data
}
\usage{
batcheffectcheck(cross, cutoff=2, minimumvariance = 0.75)
}

\arguments{
  \item{cross}{ TODO }
  \item{cutoff}{ TODO }
  \item{minimumvariance}{ minimumvariance threshold, phenotypes with lower variance are excluded }
}
\details{
  TODO
}
\value{
	A list with in each cell a batch group usable by batcheffectcorrect to correct the group means
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
		\item \code{\link{batcheffectcorrect}} - Correct batcheffects in data by normalizing group means
	}
}
\examples{
	#TODO
}

\keyword{methods}
