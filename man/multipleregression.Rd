\name{multipleregression}
\alias{multipleregression}
\title{multipleregression - multiple lineair regression on a designmatrix }
\description{
	multiple lineair regression on a designmatrix using expectation maximization and LU decomposition
}
\usage{
  multipleregression(designmatrix,y,weight=rep(1,nrow(designmatrix)),nullmodellayout=rep(1,ncol(designmatrix)),verbose=FALSE)
}

\arguments{
  \item{designmatrix}{ Designmatrix with samples in rows, and variables in columns }
  \item{y}{ Response vector (length: nrow) }
  \item{weight}{ Weight of the samples (length: nrow) }
  \item{nullmodellayout}{ which variables should be set to 0.0 to create the nullmodel (length: ncol) }
  \item{verbose}{ Be verbose, reports back model likelihood , em itterations }
}
\details{
	None Yet
}
\value{
  Likelihood of the model versus the NULL model ( Constant estimated but variable parameters set to 0.0 )
}

\references{
  TODO
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}

\seealso{
  \itemize{
    \item \code{\link{crosstocontrastlist}} - Create a contrastlist from a cross object
    \item \code{\link{backwardeliminate}} - Create a QTL model by backward elimination of cofactors
    \item \code{\link{contrastqtlmapping}} - QTL mapping routine using multiple linear regression on contrasts
    \item \code{\link{batcheffectcheck}} - Checks for batch effects in the data
		\item \code{\link{batcheffectcorrect}} - Correct batcheffects in data by normalizing group means
  }
}

\examples{
  x <- matrix(runif(100),10,1)
  y <- x[,1] + runif(2)
  o <- multipleregression(x,y)
  o
  #The example from LM
  ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
  trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
  group <- gl(2,10,20, labels=c("Ctl","Trt"))
  weight <- c(ctl, trt)
  #With intercept
  multipleregression(t(t(as.numeric(group))),y=weight)
  
  #Without intercept
  multipleregression(t(t(as.numeric(group))),y=weight,nullmodellayout=c(0))
}

\keyword{methods}
