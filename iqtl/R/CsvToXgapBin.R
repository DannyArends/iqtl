#
# CsvToXgapBin.R
#
# copyright (c) 2011, Danny Arends
# last modified dec, 2011
# first written dec, 2011
# 
# R functions: CsvToXgapBin
#

CsvToXgapBin <- function(name = "ExampleData", investigation = "ExampleInv", rowtype= "Metabolite", coltype= "Individual", valuetype="Decimal", file="metab.txt", verbose = TRUE){
  if(file.exists("CsvToBin.jar")){
    if(verbose) cat("WARNING: 'CsvToBin.jar' exists and will be re-downloaded\n")
    file.remove("CsvToBin.jar")
  }
  csvtobin <- getBinaryURL("http://www.molgenis.org/svn/standalone_tools/jars/CsvToBin.jar")
  writeBin(csvtobin,"CsvToBin.jar")
  outname <- gsub(".txt",".bin",file)
  if(file.exists(outname)){
    if(verbose) cat("WARNING: Output file '",outname,"'exists and will be overwitten\n")
    file.remove(outname)
  }
  command <- paste("java -jar CsvToBin.jar",name,investigation,rowtype,coltype,valuetype,file)
  if(verbose) cat("CMD:", command,"\n")
  res <- system(command,show.output.on.console=verbose)
  if(res != 0){
    stop("Error executing CsvToBin")
  }
}
