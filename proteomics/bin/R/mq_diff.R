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
        help="Experimntal design template", metavar="character"),
make_option(c("-f", "--file"), type="character", default=NULL,
        help="Path to the quantification file", metavar="character"),
make_option(c("-o","--output"), type="character", default=NULL,
        help="Path to the output folder", metavar="character")) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
outdir=paste(opt$o,'/',sep='')
dir.create(outdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
path=opt$f
exp_design=opt$d

source(exp_design)

data <- read.csv(path)

#ratios <- t(t(data[, cols])/colSums(data[, cols]))

intensities <- t(t(data[, cols]))

data[, cols] <- lapply(data[, cols], function(x){replace(x, x == 0,  NA)})

data[, cols] <- lapply(data[, cols], function(x){ log2(x)})

rownames(data) <- data$Identifier

Identifier <- data$Identifier
#names <- paste(Identifier, '(', data$PeptideCount, ' peptides, ', data$MS.MS.Count,' msms)',sep='')

data <- data[,cols]
data$Identifier <- Identifier
orig_data <- data
data <- data[,cols]

fpie <- function( df , names, valcols, outfile) {
    slices <- rowMeans(df[,valcols] )
    slices <- slices/sum(slices) * 100
    slices <- round(slices, 2)
    lbls <- paste(names, " ", slices, ' %', sep="")
    df$label <- t(c(lbls))
    df$slices <- t(c(slices))
    df$label[df$slices < 0.5] <- ""
    jpeg(outfile, width=1000,height=900)
    par(mar=c(6,12,6,12)+.1)
    pie(df$slices, labels = df$label,  main="Average summed intensity\nby identified taxon (%)")
    dev.off() }

out <- paste(outdir, 'mean_intensity_all.jpeg', sep='')
fpie(intensities, Identifier, cols, out)
f <- f # Defined in experimental design template
refmap <- data.frame(f, cols)

comparisons <- colnames(contrast.matrix)

groups <- rownames(contrast.matrix)

group_variance_values <- vector("list", length(groups))
group_variance_names <- vector("list", length(groups))

#variances <- data.frame()

variances <- data.frame(matrix(, nrow=nrow(orig_data), ncol=0))
print(length(orig_data))
for (group in groups){
    groupdf <- refmap[refmap$f==group,]
    groupcols <- as.numeric(rownames(groupdf))
    out <- paste(outdir, 'mean_intensity_',group,'.jpeg' , sep='')
    fpie(intensities, Identifier, groupcols, out) 
    print(groupcols)
    variances <- cbind(variances, var = apply(orig_data[,groupcols], 1, function(x) var(na.omit(x))))
    names(variances)[names(variances) == 'var'] <- group
}
dir.create(paste(outdir,'/group_variance',sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")
write.table(variances, paste(outdir,'/group_variance/group_variance.txt',sep=''), sep='\t', row.names=TRUE)

######################
# Let's make a SPLOM #
######################

print('Creating a SPLOM')
splomdata <- orig_data
splom_path=paste(outdir,'splom/',sep='')
dir.create(splom_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")

splom <- function (df, valcols, labels, path, name) {
    png(paste(path,name,sep=''),units="in",width=11,height=8.5,res=300)
    df$temp <- labels
    df <- df[df$temp!= '' , ]  
    labels <- df$temp
    df <- df[, valcols]
    prot_labels <- as.numeric(labels)
    prot_labels <- rev(rainbow_hcl(length(prot_labels)))[prot_labels]
    pairs(df, col = prot_labels, lower.panel = NULL, cex.labels=0.5, pch=19, cex = 0.01)
    dev.off() }

splom(splomdata, cols, splomdata$Identifier, splom_path, 'splom.png')

#############################################
# Diferential abundance analysis with LIMMA #
#############################################

fit2 <- lmFit(data,design)
fit2 <- contrasts.fit(fit2, contrast.matrix)
fit2 <- eBayes(fit2)

limma_dir = paste(outdir,'limma/',sep='')
dir.create(limma_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pval_lists <- vector("list", length(colnames(contrast.matrix)))
#qval_lists_0_005 <- vector("list", length(colnames(contrast.matrix)))


for ( i in seq_along(colnames(contrast.matrix))) {
    cntrst <-  colnames(contrast.matrix)[i]
    cntrst_ <- strsplit(cntrst, '-')
    Exposed <- cntrst_[[1]][1]
    Control <- cntrst_[[1]][2]
    print(cntrst)
    print('*') 
    table <- topTable(fit2,adjust="BH", coef=i, n=Inf)
    table <- merge(orig_data, table, by=0)
    table$Exposed <- Exposed
    table$Control <- Control
    table <- setDT(table, keep.rownames = TRUE)[]
    table <- table[with(table, order(P.Value)), ]
    sig_table <- table[ which(table$adj.P.Val < 0.05), ]
    
    sig_list <- sig_table$Row.names  
    pval_lists[[cntrst]] <- sig_list
    
    write.table(table, paste(limma_dir,'limma_', cntrst, '_intensity.txt',sep=''), sep='\t', row.names=FALSE)
    write.table(sig_table, paste(limma_dir,'limma_', cntrst, '_intensity_qval_0.05.txt',sep=''), sep='\t', row.names=FALSE)
    sig_list = paste(sig_list, '\n')
    write(sig_list, paste(limma_dir,'limma_', cntrst, '_intensity_qval_0.05_list.txt',sep=''))
    
    jpeg(paste(limma_dir, 'volcano_',cntrst, '_intensity.jpeg',sep=''))
    volcanoplot(fit2, coef=i)
    dev.off()
    }


num_groups = length(colnames(design))

pvals <- unlist(pval_lists, recursive = TRUE, use.names = FALSE)
pval_data <- data[rownames(data) %in% pvals, ]


#########################################################
# Hierarchical clustering and correlation of replicates #
######################################################### 
#print('Hierarchical clustering')

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



############
# Heatmaps #
############

library('dendextend')
heatmap_path=paste(outdir,'heatmaps/',sep='')
dir.create(heatmap_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
#print(pval_data)


hm <- function( df, heatmap_path, file,  main, xlab ) {
    some_col_func<-function(n)(colorspace::diverge_hcl(n,h=c(246, 40), c = 96, l = c(65, 90)))
    dend_r <- df %>% dist(method="man") %>% hclust(method="ward.D") %>% as.dendrogram %>% ladderize %>% color_branches(k=10)
    dend_c <- t(df) %>% dist(method="man") %>% hclust(method="ward.D") %>% as.dendrogram %>% ladderize%>% color_branches(k=num_groups)
    png(paste(heatmap_path, file, sep=''), units="in", width=11, height=8.5, res=300)
    gplots::heatmap.2(as.matrix(df),
    main = main,
    srtCol = 90,
    Rowv = dend_r,
    Colv = dend_c,
    trace="none",
    margins =c(10,17),
    key.xlab = xlab,
    denscol = "grey",
    density.info = "density",
    cexRow = 0.75,
    col = some_col_func)
    dev.off() }

# Heatmap of global data set
hm(data, heatmap_path, "heatmap_intensity_all.jpeg", "Log2(Intensity) - All", "log2(Intensity)")

if (length(rownames(pval_data)) >=2) {
png(paste(heatmap_path,'heatmap_significant.png',sep=''),units="in",width=11,height=8.5,res=300)
heatmap(as.matrix(pval_data), margins=c(10,17))
dev.off()}

if (length(rownames(pval_data)) >=2) {
    hm(pval_data, heatmap_path, "heatmap_significant.jpeg", "Log2(Intensity) - Adj. p-value < 0.05", "log2(Intensity)")
}

#######
# PCA #
#######
pca_path=paste(outdir,'PCA/',sep='')
dir.create(pca_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pc <- function( df, path, file, var.axes ) {
    png(paste(path,file,sep=''), units="in", width=8.5, height=8.5, res=300)
    pca <- prcomp(t(df), center=TRUE, scale=TRUE)
    g <- ggbiplot(pca, obs.scale = 1, var.scale = 1,
    groups=f, ellipse=TRUE, circle=TRUE, var.axes=var.axes)
    g <- g + scale_color_discrete(name = '')
    g <- g + theme(legend.direction = 'horizontal',legend.position = 'top')
    print(g)
    dev.off() }
pc(data, pca_path, "all_identified_pca.png", FALSE)
pc(data, pca_path, "all_identified_pca_labelled.png", TRUE)
if (length(rownames(pval_data)) >=2) {
    pc(pval_data, pca_path, "pval_0_05_identified_pca_labelled.png", TRUE)
    pc(pval_data, pca_path, "pval_0_05_identified_pca.png", FALSE)
}


# PCA significant
#pca_path=paste(outdir,'PCA/conditions/',sep='')
#dir.create(pca_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
#pc(cnd_pval_data, pca_path, "cnd_p_0_05_pca.png", TRUE)
#pc(cnd_pval_data_0_005, pca_path, "cnd_p_0_005_pca.png", TRUE)
#pc(cnd_pval_dup_data, pca_path, "cnd_p_0_05_mult_pca.png", TRUE)





#pc(data, pca_path, "all_identified_pca.png", FALSE)

#quit()

