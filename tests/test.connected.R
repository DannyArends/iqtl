#
# test.connected.R
#
# copyright (c) 2012 Danny Arends
# last modified Jan, 2012
# first written Oct, 2011
# 
require(IqtlHD)
result <- connected()
if(result==FALSE) stop("No connection to shared library: qtlHD")
