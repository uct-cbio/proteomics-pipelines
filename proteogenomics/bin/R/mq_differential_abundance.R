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
 
option_list = list(
make_option(c("-d", "--design"), type="character", default=NULL,
        help="Experimntal design template", metavar="character"),
make_option(c("-p", "--proteinGroups"), type="character", default=NULL,
        help="Path to the proteingGroups.txt file", metavar="character"),
make_option(c("-o","--output"), type="character", default=NULL,
        help="Path to the output folder", metavar="character")) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

outdir=opt$o
dir.create(outdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

path=opt$p
data <- read.csv(path)
orig_data <- data

exp_design=opt$d
source(exp_design)


data <- do.call(data.frame,lapply(data, function(x) replace(x, is.infinite(x),NA)))
data[, cols] <- lapply(data[, cols], function(x){replace(x, x == 0, NA)})
data[, cols] <- log2(data[cols])
data <- data[rowSums(is.na(data)) < 8, ]

rownames(data) <- paste(data$`X_reference.entries.mapped`)
orig_data <- data
data <- data[cols]
imputed.prob <- impute.MinProb(data)
data.missing <- data
data <- imputed.prob                        #changed data with imputed values here
fit2 <- lmFit(imputed.prob,design)
fit2 <- contrasts.fit(fit2, contrast.matrix)
fit2 <- eBayes(fit2)
limma_dir = paste(outdir,'limma/',sep='')
dir.create(limma_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pval_lists <- vector("list", length(colnames(contrast.matrix)))
pval_lists_0_005 <- vector("list", length(colnames(contrast.matrix)))

for ( i in seq_along(colnames(contrast.matrix))) {
    cntrst = colnames(contrast.matrix)[i]
    table <- topTable(fit2,adjust="BH", coef=i, n=Inf)
    table <- merge(orig_data, table, by=0)
    table <- setDT(table, keep.rownames = TRUE)[]
    pval_lists[[cntrst]] <- table[ which(table$P.Value < 0.05), ]$Row.names  # Get list of gene names with p value < 0.05
    pval_lists_0_005[[cntrst]] <- table[ which(table$P.Value < 0.005), ]$Row.names  # Get list of gene names with p value < 0.005

    table <- table[with(table, order(P.Value)), ]
    write.table(table, paste(limma_dir,'limma_', cntrst, '_iBAQ.csv',sep=''), sep=',', row.names=FALSE)
    jpeg(paste(limma_dir, 'volcano_',cntrst, '_iBAQ.jpeg',sep=''))
    volcanoplot(fit2, coef=i)
    dev.off()
    }

inter_pvals <- pval_lists[inter_comps]
inter_pvals_0_005 <- pval_lists_0_005[inter_comps]
#iso_pvals <- Reduce(intersect, iso_pvals)
inter_pvals <- unlist(inter_pvals, recursive = TRUE, use.names = FALSE)
inter_pvals_0_005  <- unlist(inter_pvals_0_005, recursive = TRUE, use.names = FALSE)
inter_pvals_duplicated <- inter_pvals[duplicated(inter_pvals)]
pval_data <- data[rownames(data) %in% inter_pvals, ]
pval_dup_data <-  data[rownames(data) %in% inter_pvals_duplicated, ]
pval_data_0_005 <- data[rownames(data) %in% inter_pvals_0_005, ]
data <- data[complete.cases(data),]
data <- na.omit(data) # listwise deletion of missing

#write.csv(data, '/researchdata/fhgfs/ptgmat003/OUTPUT/h37Rv_limma_lfq/processed_data.csv')
# plot of dataset vs dataset co-variance
#jpeg('/researchdata/fhgfs/ptgmat003/OUTPUT/h37Rv_limma_lfq/covar_LFQ.jpeg')
#pairs(data, lower.panel = NULL, cex.labels=0.5, pch=19, cex = 0.1)
#dev.off()

#
#jpeg('/researchdata/fhgfs/ptgmat003/OUTPUT/h37Rv_limma_lfq/row_dendrogram_LFQ.jpeg')
#rd <- dist(as.matrix(data), method = "man")
#rhc_fit <- hclust(rd, method = "complete")
#rdend_fit  <- as.dendrogram(rhc_fit)
#rdend  <- rotate(rdend_fit, 1:150)
#rdend <- hang.dendrogram(rdend,hang_height=0.1)
#rdend <- set(rdend, "labels_cex", 0.5)
#par(mar = c(3,3,3,7))
#plot(rdend, main = "Clustered data set", horiz =  TRUE,  nodePar = list(cex = .007))
#dev.off()

#jpeg('/researchdata/fhgfs/ptgmat003/OUTPUT/h37Rv_limma_lfq/col_dendrogram_LFQ.jpeg')
#cd <- dist(as.matrix(t(data)), method = "man")
#chc_fit <- hclust(cd, method = "complete")
#cdend_fit  <- as.dendrogram(chc_fit)
#cdend  <- rotate(cdend_fit, 1:150)
#cdend <- hang.dendrogram(cdend,hang_height=0.1)
#cdend <- set(cdend, "labels_cex", 0.5)
#par(mar = c(3,3,3,7))
#plot(cdend, main = "Clustered data set", horiz =  TRUE,  nodePar = list(cex = .007))
#dev.off()

# hierarchical clustering and correlation of replicates
reps = t(data)
hclust_path=paste(outdir,'replicate_hclust/',sep='')
dir.create(hclust_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
hclust_methods <- c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
data_dendlist <- dendlist()
for(i in seq_along(hclust_methods)) {
    jpeg(paste(hclust_path, 'replicate_cluster_', hclust_methods[i], '.jpeg', sep=''))
    tmp_dend <-  reps %>% dist(method = "man") %>% hclust(method = hclust_methods[i]) %>% as.dendrogram %>% ladderize%>% color_branches(k=3)
    labels_colors(tmp_dend) <- rainbow_hcl(6)[sort_levels_values(as.numeric(f)[order.dendrogram(tmp_dend)])]
    par(mar = c(15,3,3,3))
    plot(tmp_dend)
    dev.off()
    data_dendlist <- dendlist(data_dendlist, tmp_dend)
}

jpeg(paste(hclust_path,'replicate_cluster_corr.jpeg',sep=''))
names(data_dendlist) <- hclust_methods
cophenetic_cors <- cor.dendlist(data_dendlist)
corrplot::corrplot(cophenetic_cors, "pie", "lower")
dev.off()

#some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
some_col_func <- function(n) (colorspace::diverge_hcl(n, h = c(246, 40), c = 96, l = c(65, 90)))
heatmap_path=paste(outdir,'heatmaps/',sep='')
dir.create(heatmap_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")

#Heatmap of total data set
dend_r <- data %>% dist(method = "man") %>% hclust(method = "ward.D") %>% as.dendrogram %>% ladderize %>% color_branches(k=10)
dend_c <- t(data) %>% dist(method = "man") %>% hclust(method = "ward.D") %>% as.dendrogram %>% ladderize%>% color_branches(k=3)
jpeg(paste(heatmap_path,'heatmap_iBAQ.jpeg',sep=''))
gplots::heatmap.2(as.matrix(data),
main = "Log2(iBAQ intensity) - all",
srtCol = 45,
Rowv = dend_r,
Colv = dend_c,
trace="none",
margins =c(8,8),
key.xlab = "log2(iBAQ intensity)",
denscol = "grey",
density.info = "density",
col = some_col_func)
dev.off()

#Heatmap only of proteins with pval < 0.05 in any comparison between the two isolates.
dend_r <- pval_data %>% dist(method = "man") %>% hclust(method = "ward.D") %>% as.dendrogram %>% ladderize %>% color_branches(k=10)
dend_c <- t(pval_data) %>% dist(method = "man") %>% hclust(method = "ward.D") %>% as.dendrogram %>% ladderize%>% color_branches(k=3)
png(paste(heatmap_path,'heatmap_iBAQ_filtered_pval_0_05.png',sep=''), units="in", width=11, height=8.5, res=300)
#some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
gplots::heatmap.2(as.matrix(pval_data),
main = "Log2(iBAQ intensity) - p<0.05",
srtCol = 90,
Rowv = dend_r,
Colv = dend_c,
trace="none",
margins =c(8,10),
key.xlab = "log2(iBAQ intensity)",
denscol = "grey",
density.info = "density",
cexRow = 0.75,
col = some_col_func)
dev.off()

#Heatmap only of proteins with pval <0.05 occurring in more than 1 comparison between the two isolates.
dend_r <- pval_dup_data %>% dist(method = "man") %>% hclust(method = "ward.D") %>% as.dendrogram %>% ladderize %>% color_branches(k=10)
dend_c <- t(pval_dup_data) %>% dist(method = "man") %>% hclust(method = "ward.D") %>% as.dendrogram %>% ladderize%>% color_branches(k=3)
png(paste(heatmap_path,'heatmap_iBAQ_filtered_pval_mult_comparisons.png',sep=''), units="in", width=11, height=8.5, res=300)
#some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
gplots::heatmap.2(as.matrix(pval_dup_data),
main = "Log2(iBAQ intensity) - p<0.05 (mult. comp.)",
srtCol = 90,
Rowv = dend_r,
Colv = dend_c,
trace="none",
margins =c(8,10),
key.xlab = "log2(iBAQ intensity)",
denscol = "grey",
density.info = "density",
cexRow = 0.75,
col = some_col_func)
dev.off()

#Heatmap only of proteins with pval < 0.005 in any comparison between the two isolates.
dend_r <- pval_data_0_005 %>% dist(method = "man") %>% hclust(method = "ward.D") %>% as.dendrogram %>% ladderize %>% color_branches(k=10)
dend_c <- t(pval_data_0_005) %>% dist(method = "man") %>% hclust(method = "ward.D") %>% as.dendrogram %>% ladderize%>% color_branches(k=3)
png(paste(heatmap_path,'heatmap_iBAQ_filtered_pval_0_005.png',sep=''), units="in", width=11, height=8.5, res=300)
#some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
gplots::heatmap.2(as.matrix(pval_data_0_005),
main = "Log2(iBAQ intensity) - p<0.005",
srtCol = 90,
Rowv = dend_r,
Colv = dend_c,
trace="none",
margins =c(8,10),
key.xlab = "log2(iBAQ intensity)",
denscol = "grey",
density.info = "density",
cexRow = 0.75,
col = some_col_func)
dev.off()

# PCA
pca_path=paste(outdir,'PCA/',sep='')
dir.create(pca_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")

png(paste(pca_path,'all_data_PCA.png',sep=''), units="in", width=11, height=8.5, res=300)
pca <- prcomp(reps, center=TRUE, scale=TRUE)
#plot(pca$x[,1:2])
g <- ggbiplot(pca, obs.scale = 1, var.scale = 1,
groups = f, ellipse = TRUE, circle=TRUE, var.axes=FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
legend.position = 'top')
print(g)
dev.off()

png(paste(pca_path,'p_0_05_PCA.png',sep=''), units="in", width=11, height=8.5, res=300)
pca <- prcomp(t(pval_data), center=TRUE, scale=TRUE)
#plot(pca$x[,1:2])
library(ggbiplot)
g <- ggbiplot(pca, obs.scale = 1, var.scale = 1,
groups = f, ellipse = TRUE, circle=TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
legend.position = 'top')
print(g)
dev.off()

png(paste(pca_path,'p_0_005_PCA.png',sep=''), units="in", width=11, height=8.5, res=300)
pca <- prcomp(t(pval_data_0_005), center=TRUE, scale=TRUE)
library(ggbiplot)
g <- ggbiplot(pca, obs.scale = 1, var.scale = 1,
groups = f, ellipse = TRUE, circle=TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
legend.position = 'top')
print(g)
dev.off()

png(paste(pca_path,'p_0_05_mult_comparisons_PCA.png',sep=''), units="in", width=11, height=8.5, res=300)
pca <- prcomp(t(pval_dup_data), center=TRUE, scale=TRUE)
#plot(pca$x[,1:2])
library(ggbiplot)
g <- ggbiplot(pca, obs.scale = 1, var.scale = 1,
groups = f, ellipse = TRUE, circle=TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
legend.position = 'top')
print(g)
dev.off()

#library('VennDiagram')
#library('gplots')
#venn_path=paste(outdir, 'venn/',sep='')
#dir.create(venn_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
#png(paste(venn_path,'p_0_05_iso_venn.png',sep=''), units="in", width=11, height=8.5, res=300)
#iso_pvals <- pval_lists[iso_contrasts]
#venn.diagram(iso_pvals,
#            fill = c("red", "green", "blue"),
#            cex = 2,
#            cat.fontface = 4,
#            lty =2,
#            fontfamily =3,
#            filename = paste(venn_path,'p_0_05_iso_venn.tiff',sep=''))
#dev.off()

