\name{lodscorestoscanone}
\alias{lodscorestoscanone}
\title{lodscorestoscanone - Converts a vector/matrix of lodscores to a scanone object }
\description{
	Converts a vector/matrix of lodscores to a scanone object
}
\usage{
  lodscorestoscanone(cross,lodscores,traitnames = NULL)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{lodscores}{ vector/matrix of trait * marker to convert to scanone object }
  \item{traitnames}{ optional vector of traitnames to be used as column names for the scanone object }
}
\value{
  scanone object with the QTL results from the lodscores vector/matrix
}

\references{
  TODO
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}

\seealso{
  TODO
}

\examples{
  #TODO
}

\keyword{methods}
