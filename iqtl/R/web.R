#####################################################################
#
# web.R
#
# copyright parent_check (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Mrt, 2009
# 
# Contains: getraw - function to get an entire webpage (in lowercase) using the Rcurl package
#
######################################################################

getraw <- function(url,verbose=TRUE,logger=FALSE) {
	ourcat(logger,"GR","INFO: Trying to read url:",url,"\n",verbose=verbose)
	txt <- tolower(RCurl::getURL(url))
	ourcat(logger,"GR","Retrieved Raw:",txt,"\n",verbose=verbose)
	txt
}
