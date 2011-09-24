#library(R.lang)
#dyn.load("Dfunctions.dll")

#setwd("E:/GBIC/Arabidopsis_tiling")
#toSend <- 1:100

#.C("Test",as.integer(length(toSend)),as.integer(toSend))

getLine <- function(filename = "exp_ann.txt", linenumber = c(10,100,100000), header=TRUE, sep='\t', quoted=TRUE){
  file <- "Dfunctions.dll"
  if(!is.loaded("loadLineFromFile")) {
    dyn.load(file)
    cat(" -Loaded", file, "\n")
  }
  output <- 1:10
  OUT <- .C("loadLineFromFile",
      as.integer(getBytes(filename)),
      as.integer(linenumber),
      as.integer(header),
      getBytes(sep),
      as.integer(quoted),
      as.integer(nchar(filename)),
      as.integer(length(linenumber)),
      as.integer(output))
  OUT
}

#testdata <- getLine()
#testdata
#paste(intToChar(testdata),collapse="")

#q("no")
