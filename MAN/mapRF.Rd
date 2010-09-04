\name{mapRF}
\alias{mapRF}
\title{mapRF - QTL mapping for QTLframework using a randomforest approach }
\description{
	QTL mapping for QTLframework using a randomforest approach
}
\usage{
	mapRF(data,model)
}

\arguments{
	\item{data}{ main datastructure for QTLframework}
	\item{model}{ main modelstructure for QTLframework}
}
\details{
	Depends on the folowing R libraries:
	\itemize{
		\item{randomForest} - Available at: CRAN
	}
}
\value{

}
\references{
	Breiman, L. (2001) Random Forests \emph{Machine Learning 45} (\bold{1}) 5-32
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{

}
\seealso{
	\itemize{
		\item \code{\link{analyseQTL}} - main QTL mapping function
		\item \code{\link{mapSM}} - QTL mapping using a single marker anova model
	}
}
\examples{
	data <- testdata()
	mapRF(data)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
