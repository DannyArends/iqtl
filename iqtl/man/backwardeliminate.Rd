\name{backwardeliminate}
\alias{backwardeliminate}
\title{backwardeliminate - Create a QTL model by backward elimination of cofactors }
\description{
	Create a QTL model by backward elimination of cofactors, this model can be used with the contrastqtlmapping method or mqmscan
}
\usage{
backwardeliminate(cross,pheno.col=1,type=0,cofactors,alpha=0.01,verbose=FALSE)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{pheno.col}{ Which phenotype from the cross object are we modeling }
  \item{type}{ Which type of contrast should we use (0 Genomewide, 1 Marker condensed) }
  \item{cofactors}{ Cofactors to consider in the full model }
  \item{alpha}{ Threshold for significance }
  \item{verbose}{ Verbose output }
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
    \item \code{\link{batcheffectcheck}} - Checks for batch effects in the data
		\item \code{\link{batcheffectcorrect}} - Correct batcheffects in data by normalizing group means
    \item \code{\link{contrastqtlmapping}} - QTL mapping routine using multiple linear regression on contrasts
    \item \code{\link{contrastqtlsignificance}} - QTL permutation routine using multiple linear regression on contrasts
	}
}
\examples{
	library(iqtl)
  #Model a single traits
  data(hyper)
  hyper <- fill.geno(hyper)
  backwardeliminate(hyper,cofactors=c(3,6,9,12,17,23,45,47,78,120))
  #Model multiple traits at once
  data(multitrait)
  multi <- fill.geno(multitrait)
  backwardeliminate(multi,17:24,cofactors=c(3,6,9,12,17,23,45,47,50,63,66))
}

\keyword{methods}
