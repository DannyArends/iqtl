\name{mapSM}
\alias{mapSM}
\title{ mapSM - QTL mapping for QTLframework using a single marker anova model }
\description{
	QTL mapping for QTLframework using a single marker anova model
}
\usage{
	mapSM(data,model)
}

\arguments{
	\item{data}{ main datastructure for QTLframework}
	\item{model}{ main modelstructure for QTLframework}
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\seealso{
	\itemize{
		\item \code{\link{analyseQTL}} - main QTL mapping function
		\item \code{\link{mapRF}} - QTL mapping using randomForest algorithm
	}
}
\examples{
	data <- testdata()
	mapSM(data)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
