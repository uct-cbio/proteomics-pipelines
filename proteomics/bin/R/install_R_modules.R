#!/usr/bin/env Rscript

source("http://bioconductor.org/biocLite.R")
biocLite()

#install.packages('Rcpp', dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages('mime', dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages('optparse', dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages('devtools', dependencies=TRUE, repos='http://cran.us.r-project.org')

list.of.packages <- c("ggplot2", "Rcpp", 'mime',  'optparse', 'devtools', 'curl')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=TRUE, repos='http://cran.us.r-project.org')

biocLite('cluster',dependencies=TRUE)
biocLite('NMF',dependencies=TRUE)
biocLite('Matrix',dependencies=TRUE)
biocLite('irlba',dependencies=TRUE)
biocLite('igraph',dependencies=TRUE)
biocLite('lintr',dependencies=TRUE)
biocLite('curl',dependencies=TRUE)
biocLite('mzR',dependencies=TRUE)

biocLite('kernlab',dependencies=TRUE)
biocLite('robustbase',dependencies=TRUE)
biocLite('mclust',dependencies=TRUE)
biocLite('MLInterfaces',dependencies=TRUE)
biocLite('S4Vectors',dependencies=TRUE)
biocLite('htmltools',dependencies=TRUE)
biocLite('preprocessCore',dependencies=TRUE)
biocLite('pRoloc',dependencies=TRUE)
biocLite('MSnbase',dependencies=TRUE)
biocLite('rTANDEM',dependencies=TRUE)
biocLite('RforProteomics', dependencies=TRUE)
