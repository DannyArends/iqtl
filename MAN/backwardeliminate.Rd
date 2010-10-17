\name{backwardeliminate}
\alias{backwardeliminate}

\title{ backwardeliminate - Create a QTL model by backward elimination of cofactors }
\description{
	Create a QTL model by backward elimination of cofactors, this model can be used with the contrastqtlmapping method or mqmscan
}
\usage{
backwardeliminate(cross,pheno.col=1,cofactors,alpha=0.10)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{pheno.col}{ Which phenotype from the cross object are we modeling }
  \item{cofactors}{ Cofactors to consider in the full model }
  \item{alpha}{ Threshold for significance }
}
\details{
  TODO
}
\value{
	A list with in each cell a batch group usable by batcheffectcorrect to correct the group means
}
\references{ 
  TODO
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  TODO
}
\seealso{
	\itemize{
    \item \code{\link{highvariancephenotypes}} - Select only phenotypes that have a high variance  
		\item \code{\link{batcheffectcorrect}} - Correct batcheffects in data by normalizing group means
	}
}
\examples{
	#TODO
}

\keyword{methods}
