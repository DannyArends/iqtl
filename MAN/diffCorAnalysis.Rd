\name{diffCorAnalysis}
\alias{diffCorAnalysis}
\title{diffCorAnalysis - Calculate differential correlation (and counts) across the genome }
\description{
	Calculate differential correlation (and counts) across the genome
}
\usage{
diffCorAnalysis(cross, minimumvariance=0.01, difCorThreshold=0.4, significant = 5, method="pearson", directory="output", doplot=FALSE, writefile=FALSE, saveRdata=FALSE, verbose=TRUE)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
  \item{minimumvariance}{  Traits with less variance then minimumvariance are removed }
  \item{difCorThreshold}{ Changes in correlation >= difCorThreshold are counted in the difCntMatrix }
  \item{significant}{ Only plot probes which show differential correlation with more then significant other probes }
  \item{method}{ Method of calculating correlation }
  \item{directory}{ Directory in which to save the output }
  \item{doplot}{ Create heatmaps per marker showing differentail correlation across the genotypes }
  \item{writefile}{ Write the differential count matrix file to disk }
  \item{saveRdata}{ Save the correlation matrices as Rdata object to the directory }
  \item{verbose}{ Be verbose }
}
\details{
First thing the analysis function does is get rid of traits with no/low expression variation. Then it uses the 
\code{\link[iqtl]{differentialCorrelation}} routine to analyse all the markers it also print out the DifCor 
heatmaps in the specified directory. The counted differential expression data is saved to disk into the file 
specified by the writefile parameter
\itemize{
  \item Note: Does all the markers one by one (optimized to use 2 cores)
  \item Note: The difCor object in memory is very large. 
  \item Note: Based on the amount of traits and markers this could take a LONG time.
}
}
\value{
	A difCntMatrix object holding per row scores (which is also saved to disk) is returned invisible.
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
\examples{
	#TODO
}

\keyword{methods}
