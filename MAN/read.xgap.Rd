\name{read.xgap}
\alias{read.xgap}
\title{read.xgap - Highspeed reading of binairy data }
\description{
	Highspeed reading of binairy data structured in XGAP format from a Molgenis database or a Local file
}
\usage{
  read.xgap(location=NULL)
}

\arguments{
  \item{location}{ location (file or URL) where the binairy XGAP database file is located }
}
\details{
	Depends on the folowing R libraries: RCurl
}
\value{
  Datamatrix holding the information in the binary file
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
iqtldatadir <- paste(installed.packages()[which(rownames(installed.packages())=="iqtl"),"LibPath"],"/iqtl/data",sep="")
setwd(iqtldatadir)
r <- read.xgap("bindata.bin")
r <- read.xgap("bindataL.bin")
r <- read.xgap("bintdata.bin")
r <- read.xgap("bintdataL.bin")
}
