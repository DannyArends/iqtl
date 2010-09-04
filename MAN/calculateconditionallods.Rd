\name{calculateconditionallods}
\alias{calculateconditionallods}
\title{calculateconditionallods - Peakdetection used in inferenetwork }
\description{
	Basic peakdetection algorithm used in the inferenetwork routine
}
\usage{
  calculateconditionallods(traits, genotypes, toinfer)
}

\arguments{
  \item{traits}{ 
    Trait matrix: : cols=traits, rows=individuals
  }
  \item{genotypes}{
    Genotype matrix: cols=markers, rows=individuals
  }
  \item{toinfer}{
    List holding the locations to calculate conditional LOD scores on per trait
  }  
}
\details{
  Calculates conditional LOD scores
}

\value{
  Returns a list of matrices, holding the conditional LOD scores between the traits
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
