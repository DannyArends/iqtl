\name{find.parent}
\alias{find.parent}
\title{find.parent - Finds parents for a single individual using genetic distance measurements }
\description{
	Finds parents for a single individual in the output of the \link{score.kinship} function
}
\usage{
find.parent(res, individual, cutoff)
}

\arguments{
  \item{res}{ Result from the function score.kinship }
  \item{individual}{ Which individual are we analysing }
  \item{cutoff}{ Cutoff value (cutoff*std-dev) }
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}

\seealso{
	\itemize{
		\item \code{\link{score.kinship}} - Calculates kinshipscores based on SNP markerdata
		\item \code{\link{find.parents}} - Finds tripplets of individuals with parents from SNP markerdata
		\item \code{\link{kinship.network}} - Created a .SIF network with kinship relations between individuals
	}
}

\keyword{methods}
