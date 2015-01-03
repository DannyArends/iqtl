\name{deploy.java}
\alias{deploy.java}
\title{deploy.java - Deploy the java QTL viewer }
\description{
	Deploy the java QTL viewer
}
\usage{
	deploy.java(multiresult, cross, location)
}

\arguments{
  \item{multiresult}{ An object of class \code{mqmmulti} holding qtl results. See \code{\link[qtl]{mqmscanall}} for details. }
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details, by R/qtl }
  \item{location}{ Location to store the output files, if missing defaults to \code{\link[base]{tempdir}} }
}
\details{
  Provided upon request
}
\value{
  No return, several files are copied to location using datafile.java the data.js is created at location.
  Firefox is started when the executable is specified
}
\references{ 
  None yet
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  Only valid \code{mqmmulti} objects can be used, the objects from a mutiple QTL using scanone is not yet supported
}
\seealso{
	\itemize{
    \item \code{\link{datafile.java}} - datafile.java used internally by deploy.java
  }
}
\examples{
  require(iqtl)
  data(multitrait)
  data(locations)
  data(multitraitQTLs)
  multiloc <- addloctocross(multitrait,locations)

  # deploy.java(multitraitQTLs,cross=multiloc)
}

\keyword{iplot}
