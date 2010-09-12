\name{mqmmultitoscanone}
\alias{mqmmultitoscanone}
\title{mqmmultitoscanone - Converts a mqmmulti object to a scanone object }
\description{
	Converts a mqmmulti object to a scanone object
}
\usage{
  mqmmultitoscanone(result)
}

\arguments{
  \item{result}{ 
    Result object from \code{\link{mqmscanall}} 
  }
}
\details{
	The models inside the mqmmultiobject are not included in the resulting scanone object.
}
\value{
  scanone object with the QTL results from the mqmmultiobject
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
