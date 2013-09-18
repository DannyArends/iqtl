\name{scanonetosif}
\alias{scanonetosif}
\title{scanonetosif - scanone object to cytoscape sif file}
\description{
	Convert a scanone object to cytoscape sif file
}
\usage{
	scanonetosif(results, lodcutoff=3.0, verbose=TRUE)
}

\arguments{
  \item{results}{ An object of class \code{scanone}. See \code{\link[qtl]{scanone}} for details, by R/qtl }
  \item{lodcutoff}{ Threshold above which associations are being analyzed }
  \item{verbose}{ Display output on the conversion }
}
\details{
  Provided upon request
}
\value{
   No return, files for use with cytoscape are created in the workdir
}
\references{ 
  None yet
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  Only valid \code{cross} objects can be used
}
\seealso{
	\itemize{
    \item \code{\link{scanonetomqmmulti}} - Convert a scanone object to mqmmulti object
    \item \code{\link{qtlmatrixtomqmmulti}} - Convert a raw numeric QTL matrix to mqmmulti object
  }
}
\examples{
#TODO
}

\keyword{iplot}
