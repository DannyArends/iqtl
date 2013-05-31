#
#
# zzz.R
#
# copyright (c) 2010, Danny Arends
# last modified mrt, 2010
# first written mrt, 2010
# 
# R functions: .First.lib
#
#

.First.lib <- function(lib, pkg){
  packageStartupMessage("- Loading package Iqtl")
  library.dynam("iqtl", pkg, lib)
  .has_qtlHD <- TRUE
  tryCatch(
    library.dynam("qtlHD", pkg, lib),
    error = function(e){
     .has_qtlHD <<- FALSE
    })    
  assign(".has_qtlHD", .has_qtlHD, envir = .IqtlHDEnv)
  if(!get(".has_qtlHD", envir = .IqtlHDEnv)) packageStartupMessage(.has_qtlHD_warnmsg)
}

#Global package variables
.IqtlHDEnv <- new.env()

get_IqtlHDEnv <- function(){.IqtlHDEnv}
has_qtlHD <- function(){ get(".has_qtlHD", envir = .IqtlHDEnv) }

.has_qtlHD_warnmsg <- "- qtlHD not available, your submodules need to be initialized"

# end of zzz.R
