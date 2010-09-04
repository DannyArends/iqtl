\name{scanRF}
\alias{scanRF}
\title{ scanRF - QTL scanning using randomForest algorithm }
\description{
	Older version, please see the new mapRF function
}
\usage{
scanRF(cross = NULL, pheno.col = 1, filename = NULL, ...)
}

\arguments{
  \item{cross}{ qtl \code{cross} object }
  \item{pheno.col}{ Which phenotype (or multiple) from the cross object do we want 2 analyse }
  \item{filename}{ Save result to a file (file=<filename>) }
  \item{\dots}{ Parameters passed through to the randomForest algorithm }
}
\details{
	Depends on the folowing R libraries:
	\itemize{
		\item \link[randomForest]{randomForest} - randomForest
		\item \link[qtl]{qtl} - qtl
	}
}
\value{

}
\references{
	\itemize{
	\item Broman, K. W., Wu, H., Sen, 'S. and Churchill, G. A. (2003) R/qtl: QTL mapping in experimental crosses. \emph{Bioinformatics} (\bold{19}) 889–890
	\item Breiman, L. (2001) Random Forests \emph{Machine Learning 45} (\bold{1}) 5-32
	}
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{

}
\seealso{

}
\examples{
	require(qtl)
	data(hyper)
	result <- scanRF(cross=hyper,pheno.col=1)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
