###
# \file createRrepos.R
#
# Copyright (c) 2010, Danny Arends, p256802
# last modified May, 2011
# first written May, 2011
# 
# A copy of the GNU General Public License, version 3, is available
# at http://www.r-project.org/Licenses/GPL-3
# 
# Contains R functions: createRrepos
###
createRrepos <- function(to="./"){
  library(tools)
  setwd(to)
  if(!file.exists("GoogleFinancial_0.1.0.tar")) system("setpath r && R CMD build e:/github/Rpackages/Google-Financial")
  if(!file.exists("GoogleFinancial_0.1.0.zip")) system("setpath r && R CMD INSTALL --build e:/github/Rpackages/Google-Financial")
  
  if(!file.exists("iqtl_0.1.tar")) system("setpath r && R CMD build e:/github/Rpackages/iqtl")
  if(!file.exists("iqtl_0.1.zip")) system("setpath r && R CMD INSTALL --build e:/github/Rpackages/iqtl")
  
  if(!file.exists("pheno2geno_0.4.4.tar")) system("setpath r && R CMD build e:/github/Rpackages/phenotypes2genotypes")
  if(!file.exists("pheno2geno_0.4.4.zip")) system("setpath r && R CMD INSTALL --build e:/github/Rpackages/phenotypes2genotypes")
  
  if(!file.exists("MetaNetwork_1.0-0.tar")) system("setpath r && R CMD build e:/github/Rpackages/MetaNetwork")
  if(!file.exists("MetaNetwork_1.0-0.zip")) system("setpath r && R CMD INSTALL --build e:/github/Rpackages/MetaNetwork")
  
  if(!file.exists("qtl_1.22-2.tar")) system("setpath r && R CMD build e:/github/Rpackages/rqtl-mqm")
  if(!file.exists("qtl_1.22-2.zip")) system("setpath r && R CMD INSTALL --build e:/github/Rpackages/rqtl-mqm")
  write_PACKAGES()
}

#createRrepos()
#install.packages(contriburl="http://www.dannyarends.nl/R/")
