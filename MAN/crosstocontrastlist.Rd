\name{crosstocontrastlist}
\alias{crosstocontrastlist}
\title{crosstocontrastlist - Create a contrastlist from a cross object }
\description{
	Function to create a contrastlist from a cross object
}
\usage{
	crosstocontrastlist(cross,type=0,verbose=FALSE)
}

\arguments{
	\item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
	\item{type}{ Which type of contrast should we use (0 Genomewide, 1 Marker condensed) }
  \item{verbose}{ Verbose output }
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\seealso{
	\itemize{
    \item \code{\link{backwardeliminate}} - Create a QTL model by backward elimination of cofactors
    \item \code{\link{contrastqtlmapping}} - QTL mapping routine using multiple linear regression on contrasts
    \item \code{\link{contrastqtlsignificance}} - QTL permutation routine using multiple linear regression on contrasts
    \item \code{\link{batcheffectcheck}} - Checks for batch effects in the data
		\item \code{\link{batcheffectcorrect}} - Correct batcheffects in data by normalizing group means
	}
}
\examples{
  library(iqtl)
  data(hyper)
  hyper <- fill.geno(hyper)
  contrastlist <- crosstocontrastlist(hyper)
}
\keyword{datagen}
