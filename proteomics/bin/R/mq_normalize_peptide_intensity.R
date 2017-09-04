#!/usr/bin/env Rscript

library('limma')
library('ggbiplot')
library('gplots')
library('dendextend')
library('imputeLCMD')
library('pvclust')
library('colorspace')
library('data.table')
library('qvalue')
library("optparse")
library('MSnbase')
library("dplyr")
 
option_list = list(
make_option(c("-d", "--design"), type="character", default=NULL,
        help="Experimental design template", metavar="character"),
make_option(c("-p", "--peptides"), type="character", default=NULL,
        help="Path to the peptides.txt file", metavar="character"),
make_option(c("-o","--output"), type="character", default=NULL,
        help="Path to the output folder", metavar="character")) 

opt_parser = OptionParser(option_list=option_list);

opt = parse_args(opt_parser);

outdir=paste(opt$o,'/',sep='')

dir.create(outdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

path=opt$p

data <- read.csv(path, sep='\t')

exp_design=opt$d

source(exp_design)

data$Identifier <- data$Sequence

rownames(data) <- data$Identifier

orig_data <- data

#print(cols)

#print(data$Identifier)

data[, cols] <- lapply(data[,cols], function(x) replace(x, is.infinite(x), NA))
data[, cols] <- lapply(data[, cols], function(x){replace(x, x == 0,  NA)})
data[, cols] <- log2(data[, cols])
#data <- data[rowSums(is.na(data[,cols])) < length(cols)/2, ]
#data <- data[rowSums(is.na(data[,cols])) < 1,]

#######################################################
# Create MSnBase object, normalization and imputation #
#######################################################

print('Creating msnbase object')
msnbase_path=paste(outdir,'msnbase/',sep='')
dir.create(msnbase_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")

msnpath = paste(outdir, "msnbase/cleaned.csv",sep='')
write.csv(data, file=msnpath)

ecol <- cols
fname <- "Identifier"
eset <- readMSnSet2(msnpath, ecol, fname)
eset@phenoData$sampleNames <- cols
eset@phenoData$sampleGroups <- f

png(paste(msnbase_path,'boxplots_unnormalized.png',sep=''),units="in", width=11, height=8.5, res=300)
par(mfrow = c(2, 1))
boxplot(exprs(eset), notch=TRUE, col=(c("gold")), main="Samples", ylab="peptide log2(Intensity)", las=2) 
dev.off()

x.nrm <- normalise(eset, "quantiles")
x.imputed <- impute(x.nrm, method = "QRILC")

png(paste(msnbase_path,'boxplots_normalized.png',sep=''),units="in",width=11,height=8.5,res=300)
par(mfrow = c(2, 1))
boxplot(exprs(x.imputed), notch=TRUE, col=(c("gold")), main="Samples", ylab="peptide log2(Intensity)", las=2) 
dev.off()

png(paste(msnbase_path,'all_data_heatmap_normalized.png',sep=''),units="in",width=11,height=8.5,res=300)
heatmap(exprs(x.imputed), margins=c(10,17))
dev.off()

data <- ms2df(x.imputed)
imputedpath = paste(outdir, "msnbase/imputed.csv",sep='')
write.csv(data, file= imputedpath)

quit()
