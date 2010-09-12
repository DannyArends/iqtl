\name{deploy.svg}
\alias{deploy.svg}
\title{deploy.svg - Basic example of adding up two integers in C }
\description{
	Basic example of adding up two integers in C
}
\usage{
	deploy.svg(multiresult,cross=NULL,location="d:/",FireFox="C:/Program Files/Mozilla Firefox/firefox.exe")
}

\arguments{
  \item{multiresult}{ An object of class \code{mqmmulti} holding qtl results. See \code{\link[qtl]{mqmscanall}} for details. }
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details, by R/qtl }
  \item{location}{ location to store the output files }
  \item{FireFox}{ Location of FireFox WebBrowser }
}
\details{
  Provided upon request
}
\value{
  No return, several files are copied to location using datafile.svg the data.js is created at location.
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
    \item \code{\link{datafile.svg}} - datafile.svg used internally by deploy.svg
  }
}
\examples{
  require(iqtl)
  data(multitrait)
  data(locations)
  data(multitraitQTLs)
  multiloc <- addloctocross(multitrait,locations)
  
  deploy.svg(multitraitQTLs,cross=multiloc)
}

\keyword{iplot}
