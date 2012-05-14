#
# zzz.R
#
# copyright (c) 2012 Danny Arends
# last modified Jan, 2012
# first written Jan, 2012
# 

#Global package variables
.IqtlHDEnv <- new.env()

get_IqtlHDEnv <- function(){.IqtlHDEnv}
has_qtlHD <- function(){ get(".has_qtlHD", envir = .IqtlHDEnv) }

.has_qtlHD_warnmsg <- "- qtlHD not available, your submodules need to be initialized"

#Package loading
.onAttach <- function(lib, pkg){
  packageStartupMessage("- Loading package IqtlHD")
  .has_qtlHD <- TRUE
  tryCatch(
    library.dynam("qtlHD", pkg, lib),
    error = function(e){
     .has_qtlHD <<- FALSE
    })    
  assign(".has_qtlHD", .has_qtlHD, envir = .IqtlHDEnv)
  if(!get(".has_qtlHD", envir = .IqtlHDEnv)) packageStartupMessage(.has_qtlHD_warnmsg)
}
