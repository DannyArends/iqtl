\name{datafile.java}
\alias{datafile.java}
\title{datafile.java - mqmmulti object to datafile }
\description{
  mqmmulti object to datafile, that can be used by \code{\link{deploy.java}}
}
\usage{
	datafile.java(multiresult=NULL, cross=NULL)
}

\arguments{
  \item{multiresult}{ An object of class \code{mqmmulti} holding qtl results. See \code{\link[qtl]{mqmscanall}} for details. }
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details, by R/qtl }
}
\details{
  Provided upon request
}
\value{
  No return, a file names data.js is created in the current workdirectory. 
  The cross object is only used to extract trait genomic location from.
  When not present the CisTrans plotting option in the viewer is not available.
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
    \item \code{\link{deploy.java}} - deploy.java deploy interactive java QTL viewer
  }
}
\examples{
  require(qtl)
  data(multitrait)
  data(locations)
  data(multitraitQTLs)
  multiloc <- addloctocross(multitrait,locations)

  datafile.svg(multitraitQTLs)
}
