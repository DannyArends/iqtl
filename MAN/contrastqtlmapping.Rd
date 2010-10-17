\name{ contrastqtlmapping }
\alias{ contrastqtlmapping }
\title{ contrastqtlmapping - QTL mapping routine using multiple linear regression on contrasts }
\description{
	QTL mapping routine using multiple linear regression on contrasts
}
\usage{
	contrastqtlmapping(cross,pheno.col=1,type=0,cofactors=NULL,verbose=TRUE)
}

\arguments{
	\item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
	\item{pheno.col}{ Which phenotype from the cross object are we mapping }
	\item{type}{ Which type of contrast should we use (0 Genomewide, 1 Marker condensed) }
  \item{cofactors}{ Cofactors to consider during mapping }
  \item{verbose}{ Verbose output }
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\seealso{
	\itemize{
    \item \code{\link{batcheffectcheck}} - Checks for batch effects in the data
		\item \code{\link{batcheffectcorrect}} - Correct batcheffects in data by normalizing group means
    \item \code{\link{backwardeliminate}} - Create a QTL model by backward elimination of cofactors
	}
}
\examples{
  library(iqtl)
  data(hyper)
  hyper <- fill.geno(hyper)
  model <- backwardeliminate(hyper,cofactors=c(3,6,9,12,17,23,45,47,78,120))
  plot(contrastqtlmapping(hyper,cofactors=model))
}
\keyword{datagen}
