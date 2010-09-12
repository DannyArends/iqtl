\name{mqmmodelsasnetwork}
\alias{mqmmodelsasnetwork}
\title{mqmmodelsasnetwork - Saves a mqmmulti object as two .SIF networkfiles }
\description{
	Saves a mqmmulti object as a SIF network, this format can be visualized by using an external viewer (e.g. Cytoscape)
}
\usage{
  mqmmodelsasnetwork(cross,result)
}

\arguments{
  \item{cross}{ 
    An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{result}{
    Result object from \code{\link{mqmscanall}}
  }
}
\details{
	Two SIF files are created:
  \itemize{
    \item \code{QTLnetwork.sif} - Contains the network together with the edge attributes
    \item \code{QTLnodes.sif} - Node attribute file, containing the type of node
  }
}
\value{
	Two SIF files are created:
  \itemize{
    \item \code{QTLnetwork.sif} - Contains the network together with the edge attributes
    \item \code{QTLnodes.sif} - Node attribute file, containing the type of node
  }
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
