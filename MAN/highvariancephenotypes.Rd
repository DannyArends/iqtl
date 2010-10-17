\name{highvariancephenotypes}
\alias{highvariancephenotypes}

\title{ highvariancephenotypes - Select only phenotypes that have a high variance}
\description{
	Select only phenotypes that have a high variance
}
\usage{
highvariancephenotypes(cross, minimumvariance = 0.75)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{minimumvariance}{ minimumvariance threshold, phenotypes with lower variance are excluded }
}
\details{
  TODO
}
\value{
	Matrix phenotype values with a variance above the minimumvariance threshold
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
    \item \code{\link{batcheffectcheck}} - Checks for batch effects in the data
		\item \code{\link{batcheffectcorrect}} - Correct batcheffects in data by normalizing group means
	}
}
\examples{
	#TODO
}

\keyword{methods}
