\name{plotlodscorematrix}
\alias{plotlodscorematrix}
\title{plotlodscorematrix - Graphical representation showing QTL profile and two significance levels }
\description{
	Graphical representation showing QTL profile and two significance levels
}
\usage{
	plotlodscorematrix(cross, result, lodscorematrix)
}

\arguments{
	\item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
	\item{result}{ Result of the QTL scan }
  \item{lodscorematrix}{ matrix with each row a QTl profile obtained during permutation }
}
\details{
	Graphical representation showing QTL profile and two significance levels
}
\value{
	Plotting routine, no return value
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\seealso{
	\itemize{
    \item \code{\link{contrastqtlmapping}} - QTL mapping routine using multiple linear regression on contrasts
    \item \code{\link{contrastqtlsignificance}} - QTL permutation routine using multiple linear regression on contrasts
  }
}
\examples{
  #TODO
}

\keyword{hplot}
