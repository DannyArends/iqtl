\name{inferenetwork}
\alias{inferenetwork}
\title{inferenetwork - Causal inference of genetic networks using regression }
\description{
	Causal inference of genetic networks using regression
}
\usage{
  inferenetwork(traits, genotypes, name = "main", QTLresults = NULL,
                 lodqtl = 4, loddrop = 8)
}

\arguments{
  \item{traits}{ Trait matrix }
  \item{genotypes}{ Genotypes matrix }
  \item{name}{ Name for the results }
  \item{QTLresults}{ Pre-calculated QTL results }
  \item{lodqtl}{ Above what LOD score do we call a QTL significant }
  \item{loddrop}{ How much LOD schould be lost when regressing out T1 on T2 before calling a causal link }  
}
\details{
	Generates a SIF file holding the causaly infered network
}
\value{
  inferenceprofiles calculated during the network inference
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
data(multitrait)
res <- inferenetwork(pull.pheno(multitrait)[,1:5],pull.geno(multitrait))
}
