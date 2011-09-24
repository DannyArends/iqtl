#
#Functions inside D shared library
#First run: ./D/compile.bat for .dll
#First run: ./D/compile.sh for .so
#

StringArrayToD <- function(toSend = c("Hallo Wereld!","R->D->R")){
  file <- "Dfunctions.dll"
  dfunct <- "CharTest"
  if(!is.loaded(dfunct)){
    dyn.load(file)
    cat(" -Loading:", file, "\n")
  }
  if(is.loaded(dfunct)){
    .C("CharTest",as.integer(length(toSend)), as.integer(nchar(toSend)),as.character(toSend))
  }else{
    stop("Function",dfunct,"not in",file)
  }
}

getLine <- function(filename = "exp_ann.txt", linenumbers = c(10,100,100000), header=TRUE, sep='\t', quoted=TRUE){
  file <- "Dfunctions.dll"
  dfunct <- "loadLineFromFile"
  if(!is.loaded(dfunct)) {
    dyn.load(file)
    cat(" -Loaded", file, "\n")
  }
  if(is.loaded(dfunct)){
    output <- 0
    OUT <- .C("loadLineFromFile",
      as.character(filename),
      as.integer(linenumbers),
      as.integer(header),
      as.character(sep),
      as.integer(quoted),
      as.integer(length(filename)),
      as.integer(nchar(filename)),
      as.integer(length(linenumbers)),
      as.integer(output))
    OUT
  }else{
    stop("Function",dfunct,"not in",file)
  }
}

#StringArrayToD()
#testdata <- getLine()
#testdata
#q("no")
