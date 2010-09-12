\name{run_cluster}
\alias{run_cluster}
\title{ run_cluster - Main routine }
\description{
	run_cluster
}
\usage{
	run_cluster(name = "qtlTEST", investigation="",genotypes = "", phenotypes = "",totalitems=1, 
  njobs=1,dbpath = "", jobid=1, map="scanall",method="hk",model="normal",step=0)
}

\arguments{
  \item{name}{ Name of the resulting output matrix in molgenis }
  \item{investigation}{ Name of the investigation in the molgenis database, holding the genotypes and phenotypes }
  \item{genotypes}{ Name of genotype matrix }
  \item{phenotypes}{ Name of phenotype matrix }
  \item{totalitems}{ number of items in total }
  \item{njobs}{ number of njobs the user wants to use }
  \item{dbpath}{ URL location of the database }
  \item{jobid}{ ID of the whole job, so we can report back the status of the subjobs}
  \item{map}{ Global method used to infere QTL (scanall,cimall,mqmall) }
  \item{method}{ Method used to calculate LOD scores }
  \item{model}{ Assumed model of the underlying distribution }
  \item{step}{ Stepsize if we want to impute markers }

}
\details{
	Generates a folder
	Generates R and SH-files for use with PBS
	Submission to PBS
}
\value{
	No return value
}
\references{ 
	No refs yet
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Joeri v/d Velde
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}

\keyword{methods}
