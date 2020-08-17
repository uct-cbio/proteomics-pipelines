#!/usr/bin/env Rscript
library(matrixStats)
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
make_option(c("-q", "--quant_regex"), type="character", default=NULL,
        help="Quant regex", metavar="character"),
make_option(c("-p", "--peptides"), type="character", default=NULL,
        help="Path to the peptides.txt file", metavar="character"),
make_option(c("-o","--output"), type="character", default=NULL,
        help="Path to the output file", metavar="character"), 
make_option(c("-n","--normalize"), type="character", default=NULL,
        help="Path to the output file", metavar="character"),
make_option(c("-i","--impute"), type="character", default=NULL,
        help="Path to the output file", metavar="character")) 

opt_parser = OptionParser(option_list=option_list);

opt = parse_args(opt_parser);

outdir=paste(opt$o,'/',sep='')

dir.create(outdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

path=opt$p

data <- read.csv(path, sep='\t')

quant_regex=opt$q
norm_method=opt$n
impute_method=opt$i

#source(exp_design)

cols <- names(data)[grep(quant_regex,names(data))] 


rownames(data) <- data$Identifier

orig_data <- data


#data[, cols] <- lapply(data[,cols], function(x) {replace(x, is.infinite(x), NA)})
#data[, cols] <- t(t(data[, cols])/colSums(data[, cols]))
data[, cols] <- lapply(data[, cols], function(x){replace(x, x == 0,  NA)})

ceil <- max(data[,cols], na.rm=TRUE)
meds <- colMedians(as.matrix(data[,cols]), na.rm = TRUE)

max_med = max(meds, na.rm=TRUE)
#un.nrm[, cols] <- lapply(un.nrm[, cols], function(x){ log2(x)})
data <- data[rowSums(is.na(data[,cols])) < length(cols), ] # exclude rows where all are NA
#identifier <- data$Identifier
#data <- data[,cols]
#data$Identifier <- identifier
#print(rowSums(is.na(data[,cols])))
checkData(as.matrix(data[,cols]), verbose=TRUE)
#print(data[,cols])
#data <- data[rowSums(is.na(data[,cols])) < 2,]

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

#eset@phenoData$sampleGroups <- f

#un.nrm <- data
#un.nrm[, cols] <- lapply(un.nrm[, cols], function(x){ log2(x)})
png(paste(msnbase_path,'boxplots_unnormalized.png',sep=''),units="in", width=11, height=8.5, res=300)
par(mfrow = c(2, 1))
par(cex.lab=0.5) # is for y-axis
par(cex.axis=0.5) # is for x-axis

boxplot(log2(exprs(eset)), notch=TRUE, outline=FALSE, col=(c("gold")), main="Samples", ylab="log2(Intensity)", las=2) 
dev.off()

#x.nrm <- eset
print(paste("Normalizing data using: ", norm_method,sep=''))

checkData(as.matrix(exprs(eset)), verbose=TRUE)
baseline = min(exprs(eset), na.rm=TRUE)
if (norm_method != 'none') {
    x.nrm <- normalise(eset, norm_method) 
} else{
    x.nrm <- eset
}

# https://pubs.acs.org/doi/pdf/10.1021/acs.jproteome.5b00981
checkData(as.matrix(exprs(x.nrm)), verbose=TRUE)

x.imputed <- impute(x.nrm, method = impute_method, colmax=90)

minval = min(exprs(x.imputed), na.rm=TRUE)
delta <- baseline - minval
exprs(x.imputed) <- exprs(x.imputed) + delta # ensures the minimum value is the same after imputation - BPCA can yield negative values 
minval = min(exprs(x.nrm), na.rm=TRUE)
delta <- baseline - minval
exprs(x.nrm) <- exprs(x.nrm) + delta # ensures the minimum value is the same after imputation - BPCA can yield negative values 
#exprs(x.nrm) <- exprs(x.nrm) 
#x.imputed <- impute(x.nrm, method = impute_method)

#x.imputed <- x.nrm

png(paste(msnbase_path,'boxplots_normalized.png',sep=''),units="in",width=11,height=8.5,res=300)
par(mfrow = c(2, 1))
par(cex.lab=0.5) # is for y-axis
par(cex.axis=0.5) # is for x-axis
boxplot(log2(exprs(x.nrm) ) , notch=TRUE, col=(c("gold")), outline=FALSE, main="Samples", ylab="log2(Intensity)", las=2) 
dev.off()

png(paste(msnbase_path,'boxplots_imputed.png',sep=''),units="in",width=11,height=8.5,res=300)
par(mfrow = c(2, 1))
par(cex.lab=0.5) # is for y-axis
par(cex.axis=0.5) # is for x-axis
boxplot(log2(exprs(x.imputed)) , notch=TRUE, col=(c("gold")), main="Samples", outline=FALSE, ylab="log2(Intensity)", las=2) 
dev.off()

x.nrm <- x.imputed

png(paste(msnbase_path,'all_data_heatmap_normalized.png',sep=''),units="in",width=11,height=8.5,res=300)
heatmap(log2(exprs(x.nrm)), margins=c(10,17))
dev.off()

data <- ms2df(x.nrm)
#data[, cols] <- lapply(data[, cols], function(x){2^x})
#data[, cols] <- lapply(data[, cols], function(x){replace(x, x == NA,  0)})

imputedpath = paste(outdir, "msnbase/normalized.csv",sep='')
write.csv(data, file= imputedpath)

# shuffle it for permutation testing
data <- data[sample(nrow(data)),]
imputedpath = paste(outdir, "msnbase/normalized.csv.shuffled.csv",sep='')
write.csv(data, file= imputedpath)
