\name{correlation}
\alias{correlation}
\alias{matrix_correlation}
\title{correlation - Compute correlations (D programming language) }
\description{
  Compute correlations (D programming language)
}
\usage{
correlation(x, y, verbose = FALSE)
}

\arguments{
  \item{x}{ Numeric vector/matrix to calculate the correlation of }
  \item{y}{ Numeric vector to calculate the correlation between x and y. }
  \item{verbose}{ Be verbose. }
}
\details{
  This function uses a call into a shared object produced by using the D programming language.
  When it is unable to build the shared object, it uses the \code{\link[stats]{cor}} routine 
  as a fall back.
}
\value{
  Correlation vector/matrix depending on the input
}
\references{ 
  TODO
}
\author{ 
  Danny Arends \email{Danny.Arends@gmail.com}\cr
  Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  TODO
}
\examples{
  require(IqtlHD)
  matr <- matrix(runif(100*300),100,300)
  #res1 <- correlation(matr)
  #res2 <- correlation(matr[,1],matr[,2])
}

\keyword{methods}
