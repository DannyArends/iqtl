\name{getlodthreshold}
\alias{getlodthreshold}
\title{getlodthreshold - Calculate lod cut-off for a user specified significance level }
\description{
	Calculate lod cut-off for a user specified significance level out of a matrix of lodscores (from permutation)
}
\usage{
  getlodthreshold(lodscorematrix,percentage = 5)
}

\arguments{
  \item{lodscorematrix}{ matrix with each row a QTl profile obtained during permutation }
  \item{percentage}{ user specified significance level }
}
\value{
  lod cut-off for the user specified significance level
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
