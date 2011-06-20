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
createRrepos <- function(to="E:/eclipse/Production/3Dto2DApplet/websites/homepage/Rpackages"){
  library(tools)
  setwd(to)
  if(!file.exists("GoogleFinancial_0.1.0.tar")) system("setpath r && R CMD build e:/github/Google-Financial")
  if(!file.exists("GoogleFinancial_0.1.0.zip")) system("setpath r && R CMD build e:/github/Google-Financial --binary")
  
  if(!file.exists("iqtl_0.1.tar")) system("setpath r && R CMD build e:/iqtl")
  if(!file.exists("iqtl_0.1.zip")) system("setpath r && R CMD build e:/iqtl --binary")
  
  if(!file.exists("pheno2geno_0.4.4.tar")) system("setpath r && R CMD build e:/github/phenotypes2genotypes")
  if(!file.exists("pheno2geno_0.4.4.zip")) system("setpath r && R CMD build e:/github/phenotypes2genotypes --binary")
  
  if(!file.exists("MetaNetwork_1.0-0.tar")) system("setpath r && R CMD build e:/MetaNetwork")
  if(!file.exists("MetaNetwork_1.0-0.zip")) system("setpath r && R CMD build e:/MetaNetwork --binary")
  
  if(!file.exists("qtl_1.22-2.tar")) system("setpath r && R CMD build e:/github/rqtl-mqm")
  if(!file.exists("qtl_1.22-2.zip")) system("setpath r && R CMD build e:/github/rqtl-mqm --binary")
  write_PACKAGES()
}

createRrepos()
#install.packages(contriburl="http://www.dannyarends.nl/Rpackages/")