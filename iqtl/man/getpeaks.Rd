\name{getpeaks}
\alias{getpeaks}
\title{getpeaks - Peakdetection used in inferenetwork }
\description{
	Basic peakdetection algorithm used in the inferenetwork routine
}
\usage{
  getpeaks(qtlprofiles, cutoff = 4.0)
}

\arguments{
  \item{qtlprofiles}{ 
    Matrix of QTL profiles created here: cols=markers, rows=traits
  }
  \item{cutoff}{
    Above which LOD score do we call a QTL significant
  }
}
\details{
  Basic peak detection, to select the top marker in a region of markers above a certain threshold
}

\value{
  Returns a list of vectors, where the numbers inside the vectors are marker numbers found to be significant (above the cutoff).
  Each list item is a trait.
}

\references{
  TODO
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}

\seealso{
  TODO
}

\examples{
#TODO
}
