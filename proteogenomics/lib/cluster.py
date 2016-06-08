#!/usr/bin/python

import numpy as np
import scipy as sp
import pandas as pd
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import tempfile
import os
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import string
import rpy2
from rpy2 import robjects

input ='''
library(gplots)
library(dendextend)
# This is the input part
# Filtering and restructuring
all_cluster_data <- do.call(data.frame,lapply(all_cluster_data, function(x) replace(x, is.infinite(x),NA)))
filtered_all_cluster_data <- all_cluster_data[complete.cases(all_cluster_data),]
gene_id <- filtered_all_cluster_data[,2]
# Removing execess columns (Check to make sure correct columns are removed)
filtered_all_cluster_data<- filtered_all_cluster_data[,-2]
filtered_all_cluster_data<- filtered_all_cluster_data[,-1]
#filtered_all_cluster_data<- filtered_all_cluster_data[,]
drops <- c("X")
filtered_all_cluster_data <- filtered_all_cluster_data[ , !(names(filtered_all_cluster_data) %in% drops)]
print(names(filtered_all_cluster_data))
# Taking only the transcriptome data and creating a new dataset (Also changes if columns move about)
#trans_only_set <- filtered_all_cluster_data[,c(1,2,3,4)]
'''

covariance ='''
# plot of dataset vs dataset co-variance
#pairs(filtered_all_cluster_data, lower.panel = NULL, cex.labels=0.5, pch=19, cex = 0.1)
pairs(filtered_all_cluster_data, lower.panel = NULL, cex.labels=0.05, pch=19, cex = 0.1)
'''

testing = '''
# Random other code that was used in testing
transpose_set <- t(filtered_all_cluster_data)
subset <- trans_set[-c(1,2), ]
head(subset, n=10)
distance_matrix <- dist(subset, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
aclust <- hclust(distance_matrix, method = "complete")
plot(aclust)
'''


'''
# dendrogram and clustering of only transcript data
dist_transcript <- dist(trans_only_set)
hc_transcript <- hclust(dist_transcript, method = "complete")
dend_transcrip <- as.dendrogram(hc_transcript)
dend_transcrip <- rotate(dend_transcrip, 1:150)
dend_transcrip <- hang.dendrogram(dend_transcrip,hang_height=0.1)
dend_transcrip <- set(dend_transcrip, "labels_cex", 0.5)
par(mar = c(3,3,3,7))
plot(dend_transcrip, main = "Clustered transcriptome data set", horiz =  TRUE,  nodePar = list(cex = .007))
some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
gplots::heatmap.2(as.matrix(trans_only_set), 
                  main = "Heatmap for transcriptome data set",
                  srtCol = 20,
                  dendrogram = "row",
                  Rowv = dend_transcrip,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(5,0.1),      
                  key.xlab = "Log2FC",
                  denscol = "grey",
                  density.info = "density",        
                  col = some_col_func
)
'''

def cluster(data, outdir,name, xlab):
    name = name.split()
    name = '_'.join(name)
    rdf =pandas2ri.py2ri(data)
    ro.globalenv['all_cluster_data'] = rdf
    ro.r(input)
    ro.r("jpeg('{}/{}_cov.jpeg')".format(outdir, name))
    ro.r(covariance)
    ro.r('dev.off()')
    ro.r("jpeg('{}/{}_heatmap.jpeg')".format(outdir, name))

    dendrogram='''
    # dendrogram and clustering of all data
    d_filtered_all_cluster_data <- dist(filtered_all_cluster_data)
    hc_filtered_all_cluster_data <- hclust(d_filtered_all_cluster_data, method = "complete")
    dend_filtered_all_cluster_data <- as.dendrogram(hc_filtered_all_cluster_data)
    dend_filtered_all_cluster_data <- rotate(dend_filtered_all_cluster_data, 1:150)
    dend_filtered_all_cluster_data <- hang.dendrogram(dend_filtered_all_cluster_data,hang_height=0.1)
    dend_filtered_all_cluster_data <- set(dend_filtered_all_cluster_data, "labels_cex", 0.5)
    par(mar = c(3,3,3,7))
    plot(dend_filtered_all_cluster_data, main = "Clustered data set", horiz =  TRUE,  nodePar = list(cex = .007))
    some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
    head(filtered_all_cluster_data)
    gplots::heatmap.2(as.matrix(filtered_all_cluster_data),
                      main = "Heatmap for total data set",
                      srtCol = 20,
                      dendrogram = "row",
                      Rowv = dend_filtered_all_cluster_data,
                      Colv = "NA", # this to make sure the columns are not ordered
                      trace="none",
                      margins =c(5,0.1),
                      key.xlab = "{}",
                      denscol = "grey",
                      density.info = "density",
                      col = some_col_func
    )
    '''.format(xlab)
    ro.r(dendrogram)
    ro.r('dev.off()')
