\name{multipleregression}
\alias{multipleregression}
\title{multipleregression - multiple lineair regression on a designmatrix }
\description{
	multiple lineair regression on a designmatrix using expectation maximization and LU decomposition
}
\usage{
  multipleregression(designmatrix,weight,y,verbose=FALSE)
}

\arguments{
  \item{designmatrix}{ Designmatrix with samples in rows, and variables in columns }
  \item{weight}{ Weight of the samples (length: nrow) }
  \item{y}{ Response vector (length: nrow) }
  \item{verbose}{ Be verbose, reports back model likelihood , em itterations }
}
\details{
	None Yet
}
\value{
  Likelihood of the model versus the NULL model ( Constant estimated but variable parameters set to 0.0 )
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
  x <- matrix(runif(100),10,1)
  w <- rep(1,10)
  y <- x[,1] + runif(2)
  o <- multipleregression(x,w,y)
  o
}

\keyword{IO}
\keyword{database}
