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

#exit()

option_list = list(
make_option(c("-d", "--design"), type="character", default=NULL,
        help="Experimntal design template", metavar="character"),
make_option(c("-p", "--proteinGroups"), type="character", default=NULL,
        help="Path to the proteingGroups.txt file", metavar="character"),
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

rownames(data) <- data$Identifier

orig_data <- data

#data[, cols] <- lapply(data[,cols], function(x) replace(x, is.infinite(x),NA))

#data[, cols] <- lapply(data[, cols], function(x){replace(x, x == 0, NA)})

#data[, cols] <- log2(data[, cols])

#exclude missing samples 
#data <- data[rowSums(is.na(data[,cols])) < length(cols)/2, ]


#######################################################
# Create MSnBase object, normalization and imputation #
#######################################################
#print('Creating msnbase object')
#msnbase_path=paste(outdir,'msnbase/',sep='')
#dir.create(msnbase_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")

#msnpath = paste(outdir, "msnbase/cleaned.csv",sep='')
#write.csv(data, file=msnpath)

#ecol <- cols
#fname <- "Identifier"
#eset <- readMSnSet2(msnpath, ecol, fname)
#eset@phenoData$sampleNames <- cols
#eset@phenoData$sampleGroups <- f

#png(paste(msnbase_path,'boxplots_unnormalized.png',sep=''),units="in", width=11, height=8.5, res=300)
#par(mfrow = c(2, 1))
#boxplot(exprs(eset), notch=TRUE, col=(c("gold")), main="Samples", ylab="protein log2(iBAQ intensity)", las=2) 
#dev.off()

#x.nrm <- normalise(eset, "quantiles")
#x.imputed <- impute(x.nrm, method = "bpca")

#png(paste(msnbase_path,'boxplots_normalized.png',sep=''),units="in",width=11,height=8.5,res=300)
#par(mfrow = c(2, 1))
#boxplot(exprs(x.imputed), notch=TRUE, col=(c("gold")), main="Samples", ylab="protein log2(iBAQ intensity)", las=2) 
#dev.off()

#png(paste(msnbase_path,'all_data_heatmap_normalized.png',sep=''),units="in",width=11,height=8.5,res=300)
#heatmap(exprs(x.imputed), margins=c(10,17))
#dev.off()

#data <- ms2df(x.imputed)

orig_data <- data[,cols]

identifier <- data$Identifier

orig_data$Identifier <- identifier 

refids <- data$'Gene_OrderedLocusName'

orig_data$'Gene_OrderedLocusName' <- refids
orig_data$'UniProtKB_ID' <- data$"UniProtKB.ID"

orf_ids <- data$"ORF.ids...all.strains" 
orig_data$"ORF.ids...all.strains" <- orf_ids

pg_ids  <- data$'Identifier'

orig_data$'Identifier' <- pg_ids

GI <- data$"GeneID" 
orig_data$"GeneID" <- GI

orig_data$'iBAQMean' <- rowMeans( data[,cols] )

print('Removing isoform pg based on combined orf ids')
orig_data <- orig_data[order(-orig_data$iBAQMean),] # Order by mean ibaq, descending

orig_data <- orig_data[!duplicated(orig_data$"ORF.ids...all.strains" ), ]

#imputedpath = paste(outdir, "msnbase/imputed.csv",sep='')
#write.csv(orig_data, file= imputedpath)


data <- data[,cols]

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

splom(splomdata, cols, splomdata$Identifier, splom_path, 'all_proteins_splom.png')

#splom(splomdata, cols, splomdata$Protein.families, splom_path,'all_proteins_families_splom.png')
#splom(splomdata, cols, splomdata$Pathway, splom_path, 'all_proteins_pathways_splom.png')


#############################################
# Diferential abundance analysis with LIMMA #
#############################################

fit2 <- lmFit(data,design)
fit2 <- contrasts.fit(fit2, contrast.matrix)
fit2 <- eBayes(fit2)

limma_dir = paste(outdir,'limma/',sep='')
dir.create(limma_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pval_lists <- vector("list", length(colnames(contrast.matrix)))
qval_lists <- vector("list", length(colnames(contrast.matrix)))

pval_lists_0_005 <- vector("list", length(colnames(contrast.matrix)))

for ( i in seq_along(colnames(contrast.matrix))) {
    cntrst = colnames(contrast.matrix)[i]
    cntrst_ <- strsplit(cntrst, '-')
    Exposed <- cntrst_[[1]][1]
    Control <- cntrst_[[1]][2]
    #print(Exposed)
    #print(Control)
    #print('*') 
    table <- topTable(fit2,adjust="BH", coef=i, n=Inf)
    table <- merge(orig_data, table, by=0)
    
    
    table$Exposed <- Exposed
    table$Control <- Control
    
    table <- setDT(table, keep.rownames = TRUE)[]
    #print(head(table)) 

    pval_lists[[cntrst]] <- table[ which(table$P.Value < 0.05), ]$Row.names  
    qval_lists[[cntrst]] <- table[ which(table$adj.P.Val < 0.05), ]$Row.names  
    print(cntrst)
    print(qval_lists) 
    print("***")
    pval_lists_0_005[[cntrst]] <- table[ which(table$P.Value < 0.005), ]$Row.names  
    table <- table[with(table, order(P.Value)), ]
    write.table(table, paste(limma_dir,'limma_', cntrst, '_iBAQ.csv',sep=''), sep='\t', row.names=FALSE)
    jpeg(paste(limma_dir, 'volcano_',cntrst, '_iBAQ.jpeg',sep=''))
    volcanoplot(fit2, coef=i)
    dev.off()
    }

num_groups = length(colnames(design))

# Get the proteins that are differentially expressed between strains
inter_pvals <- pval_lists[inter_comps]
inter_qvals <- qval_lists[inter_comps]

inter_pvals_0_005 <- pval_lists_0_005[inter_comps]

inter_pvals <- unlist(inter_pvals, recursive = TRUE, use.names = FALSE)

print(inter_qvals)
inter_qvals <- unlist(inter_qvals, recursive = TRUE, use.names = FALSE)

print(inter_qvals)

inter_pvals_0_005  <- unlist(inter_pvals_0_005, recursive = TRUE, use.names = FALSE)
inter_pvals_duplicated <- inter_pvals[duplicated(inter_pvals)]

str_pval_data <- data[rownames(data) %in% inter_pvals, ]
str_qval_data <- data[rownames(data) %in% inter_qvals, ]
str_pval_dup_data <-  data[rownames(data) %in% inter_pvals_duplicated, ]
str_pval_data_0_005 <- data[rownames(data) %in% inter_pvals_0_005, ]

# Get the proteins that are differentially expressed between conditions
intra_pvals <- pval_lists[intra_comps]
intra_qvals <- qval_lists[intra_comps]


intra_pvals_0_005 <- pval_lists_0_005[intra_comps]

intra_pvals <- unlist(intra_pvals, recursive = TRUE, use.names = FALSE)
intra_qvals <- unlist(intra_qvals, recursive = TRUE, use.names = FALSE)

intra_pvals_0_005  <- unlist(intra_pvals_0_005, recursive = TRUE, use.names = FALSE)
intra_pvals_duplicated <- intra_pvals[duplicated(intra_pvals)]

cnd_pval_data <- data[rownames(data) %in% intra_pvals, ]
cnd_qval_data <- data[rownames(data) %in% intra_qvals, ]
cnd_pval_dup_data <-  data[rownames(data) %in% intra_pvals_duplicated, ]
cnd_pval_data_0_005 <- data[rownames(data) %in% intra_pvals_0_005, ]

data <- data[complete.cases(data),]
data <- na.omit(data) # listwise deletion of missing

#########################################################
# Hierarchical clustering and correlation of replicates #
######################################################### 
print('Hierarchical clustering')

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

jpeg(paste(hclust_path,'replicate_cluster_corr.jpeg',sep=''))
names(data_dendlist) <- hclust_methods
cophenetic_cors <- cor.dendlist(data_dendlist)
corrplot::corrplot(cophenetic_cors, "pie", "lower")
dev.off()

############
# Heatmaps #
############

library('dendextend')
#dend_r <- df %>% dist(method="man") %>% hclust(method="ward.D") %>% as.dendrogram %>% ladderize %>% color_branches(k=10)
dend_c <- t(data) %>% dist(method="man") %>% hclust(method="ward.D") %>% as.dendrogram %>% ladderize%>% color_branches(k=num_groups)
dend_r <- data %>% dist(method="man") %>% hclust(method="ward.D") %>% as.dendrogram %>% ladderize %>% color_branches(k=10)

heatmap_path=paste(outdir,'heatmaps/',sep='')

dir.create(heatmap_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")

hm <- function( df, heatmap_path, file,  main, xlab, dend_c, dend_r ) {
    some_col_func<-function(n)(colorspace::diverge_hcl(n,h=c(246, 40), c = 96, l = c(65, 90)))
    #dend_r <- df %>% dist(method="man") %>% hclust(method="ward.D") %>% as.dendrogram %>% ladderize %>% color_branches(k=10)
    #dend_c <- t(df) %>% dist(method="man") %>% hclust(method="ward.D") %>% as.dendrogram %>% ladderize%>% color_branches(k=num_groups)
    png(paste(heatmap_path, file, sep=''), units="in", width=11, height=8.5, res=300)
    gplots::heatmap.2(as.matrix(df),
    main = main,
    srtCol = 90,
    #Rowv = dend_r,
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
hm(data, heatmap_path, "heatmap_iBAQ.jpeg", "Log2(iBAQ intensity) - all", "log2(iBAQ intensity", dend_c, dend_r)
# Strain comparisons
heatmap_path=paste(outdir,'heatmaps/strains/',sep='')
dir.create(heatmap_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
title="Log2(iBAQ intensity) - strain comparisons (p < 0.05)"
hm(str_pval_data, heatmap_path, "heatmap_strains_05_iBAQ.jpeg", title, "log2(iBAQ intensity", dend_c, dend_r)
title="Log2(iBAQ intensity) - str. comps. mult. (p < 0.05)"
hm(str_pval_dup_data, heatmap_path, "heatmap_strains_05_mult_iBAQ.jpeg", title, "log2(iBAQ intensity", dend_c, dend_r)
title="Log2(iBAQ intensity) - str. comps. (p < 0.005)"
hm(str_pval_data_0_005, heatmap_path, "heatmap_strains_005_iBAQ.jpeg", title, "log2(iBAQ intensity", dend_c, dend_r)

title="Log2(iBAQ intensity) - strain comparisons (q < 0.05)"
hm(str_qval_data, heatmap_path, "heatmap_strains_q05_iBAQ.jpeg", title, "log2(iBAQ intensity", dend_c, dend_r)

# Condition comparisons
heatmap_path=paste(outdir,'heatmaps/conditions/',sep='')
dir.create(heatmap_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
title="Log2(iBAQ intensity) - condition comparisons (p < 0.05)"
hm(cnd_pval_data, heatmap_path, "heatmap_cnd_05_iBAQ.png", title, "log2(iBAQ intensity", dend_c, dend_r)
title="Log2(iBAQ intensity) - cnd. comps. mult. (p < 0.05)"
hm(cnd_pval_dup_data, heatmap_path, "heatmap_cnd_05_mult_iBAQ.png", title, "log2(iBAQ intensity", dend_c, dend_r)
title="Log2(iBAQ intensity) - cnd. comps. (p < 0.005)"
hm(cnd_pval_data_0_005, heatmap_path, "heatmap_cnd_005_iBAQ.jpeg", title, "log2(iBAQ intensity", dend_c, dend_r)

title="Log2(iBAQ intensity) - condition comparisons (q < 0.05)"
hm(cnd_qval_data, heatmap_path, "heatmap_conditions_q05_iBAQ.jpeg", title, "log2(iBAQ intensity", dend_c, dend_r)

quit()
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
    g <- g + theme(legend.direction = 'horizontal',
    legend.position = 'top')
    #print(g)
    dev.off() }

pc(data, pca_path, "all_identified_pca.png", FALSE)

# PCA strains
pca_path=paste(outdir,'PCA/strains/',sep='')
dir.create(pca_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pc(str_pval_data, pca_path, "str_p_0_05_pca.png", TRUE)
pc(str_pval_data_0_005, pca_path, "str_p_0_005_pca.png", TRUE)
pc(str_pval_dup_data, pca_path, "str_p_0_05_mult_pca.png", TRUE)

# PCA conditions
pca_path=paste(outdir,'PCA/conditions/',sep='')
dir.create(pca_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pc(cnd_pval_data, pca_path, "cnd_p_0_05_pca.png", TRUE)
pc(cnd_pval_data_0_005, pca_path, "cnd_p_0_005_pca.png", TRUE)
pc(cnd_pval_dup_data, pca_path, "cnd_p_0_05_mult_pca.png", TRUE)

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

