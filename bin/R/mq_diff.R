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
data <- data[with(data, order(-rowSums(data[,cols, drop=FALSE]))), ]

#ratios <- t(t(data[, cols])/colSums(data[, cols]))

intensities <- t(t(data[, cols]))

data[, cols] <- lapply(data[, cols], function(x){replace(x, x == 0,  NA)})

data[, cols] <- lapply(data[, cols], function(x){ log2(x)})

rownames(data) <- data$Identifier

Identifier <- data$Identifier
#names <- paste(Identifier, '(', data$PeptideCount, ' peptides, ', data$MS.MS.Count,' msms)',sep='')

#data <- data[,cols]
data$Identifier <- Identifier
orig_data <- data
data <- data[,cols]

fpie <- function( df , names, valcols, outfile) {
    df$slices <- rowSums(df[,valcols, drop=FALSE] )
    df$slices <- df$slices/sum(df$slices) * 100
    df$slices <- round(df$slices, 2)
    df$label <- paste(names, " ", df$slices, ' %', sep="")
    df$label[df$slices < 0.5] <- ""
    jpeg(outfile, width=1000,height=900)
    par(mar=c(6,12,6,12)+.1)
    pie(df$slices, labels = df$label,  main="Summed intensity (%)")
    dev.off() }

out <- paste(outdir, 'summed_intensity_all.jpeg', sep='')
fpie(intensities, Identifier, cols, out)
f <- f # Defined in experimental design template
refmap <- data.frame(f, cols)
comparisons <- colnames(contrast.matrix)
groups <- rownames(contrast.matrix)
group_variance_values <- vector("list", length(groups))
group_variance_names <- vector("list", length(groups))
#variances <- data.frame()
vardata <- data
vardata$Identifier <- Identifier
vardata$sums  <- rowSums(intensities[,cols, drop=FALSE] ) # sum the raw intensity to exclude proteins that contribute less than 1 percent
total <- sum(vardata$sums)
vardata$percentage <- lapply(vardata$sums, function(x){(x/total)*100})
vardata <- vardata[vardata$percentage >= 1,]

variances <- data.frame(matrix(, nrow=nrow(vardata), ncol=0))
dir.create(paste(outdir,'/group_variance',sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")
for (group in groups){
    groupdf <- refmap[refmap$f==group,]
    groupcols <- as.numeric(rownames(groupdf))
    out <- paste(outdir, 'summed_intensity_',group,'.jpeg' , sep='')
    fpie(intensities, Identifier, groupcols, out) 
    variances <- cbind(variances, var = apply(vardata[,groupcols, drop=FALSE], 1, function(x) var(na.omit(x))))
    names(variances)[names(variances) == 'var'] <- group 
}

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
    table <- topTable(fit2,adjust="BH", coef=i, n=Inf)
    #table <- merge(orig_data, table, by=0)
    table$Exposed <- Exposed
    table$Control <- Control
    table <- setDT(table, keep.rownames = TRUE)[]
    table <- table[with(table, order(P.Value)), ]
    sig_table <- table[ which(table$adj.P.Val < 0.05), ]
    
    sig_list <- sig_table$Row.names  
    print(cntrst)
    print(sig_list)
    print('*') 
    pval_lists[[cntrst]] <- sig_list
    
    write.table(table, paste(limma_dir,'limma_', cntrst, '_intensity.csv',sep=''), sep=',', row.names=FALSE)
    write.table(sig_table, paste(limma_dir,'limma_', cntrst, '_intensity_qval_0.05.csv',sep=''), sep=',', row.names=FALSE)
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


#dend_c <- t(data) %>% dist(method="man") %>% hclust(method="ward.D") %>% as.dendrogram %>% ladderize%>% color_branches(k=num_groups)

hm <- function( df, heatmap_path, file,  main, xlab , factors) {
    num_groups <- length(unique(factors))
    some_col_func<-function(n)(colorspace::diverge_hcl(n,h=c(246, 40), c = 96, l = c(65, 90)))
    col_colors<-rainbow_hcl(length(unique(factors)))[sort_levels_values(as.numeric(factors))]
    #dend_r <- df %>% dist(method="man") %>% hclust(method="ward.D") %>% as.dendrogram %>% ladderize %>% color_branches(k=10)
    dend_c <- t(df) %>% dist(method="man") %>% hclust(method="ward.D") %>% as.dendrogram %>% ladderize%>% color_branches(k=num_groups)
    png(paste(heatmap_path, file, sep=''), units="in", width=11, height=8.5, res=300)
    gplots::heatmap.2(as.matrix(df),
    main = main,
    srtCol = 90,
    #Rowv = dend_r,
    Colv = dend_c,
    ColSideColors = col_colors,
    trace="none",
    margins =c(15,15),
    key.xlab = xlab,
    denscol = "grey",
    density.info = "density",
    cexRow = 0.75,
    col = some_col_func,
    colCol= col_colors )
    par(lend = 1)           # square line ends for the color legend
    legend("topright",      # location of the legend on the heatmap plot
               legend = unique(factors), # category labels
               col =unique(col_colors),  # color key
                  lty= 1,             # line style
                       lwd = 10            # line width
                   )
    dev.off() }

# Heatmap of global data set
hm(data, heatmap_path, "heatmap_intensity_all.png", "Log2(Intensity) - All", "log2(Intensity)", f)

if (length(rownames(pval_data)) >=2) {
    hm(pval_data, heatmap_path, "heatmap_significant.png", "Log2(Intensity) - Adj. p-value < 0.05", "log2(Intensity)", f)
}
for ( i in seq_along(colnames(contrast.matrix))) {
    cntrst <-  colnames(contrast.matrix)[i]
    cntrst_ <- strsplit(cntrst, '-')
    fields <- c()
    levels <- c()

    for (group in cntrst_[[1]] ) {
        groupdf <- refmap[refmap$f==group,]
        groupcols <- as.numeric(rownames(groupdf))
        grouplen <- length(groupcols)
        glevels   <- rep(group, grouplen) 
        fields <- c(fields, groupcols)
        levels <- c(levels, glevels)
    }

    compdata <- data[,fields]
    sig_list <- pval_lists[[cntrst]]
    group_pval_data <- compdata[rownames(compdata) %in% sig_list, ]
    newf <- droplevels(f[fields])
    #levels(newf) <- levels
    if (length(rownames(group_pval_data)) >=2) {
        hm(group_pval_data, heatmap_path, paste(cntrst, "heatmap_significant.png",sep='_'), "Log2(Intensity) - Adj. p-value < 0.05", "log2(Intensity)", newf)
    }
}


#######
# PCA #
#######
pca_path=paste(outdir,'PCA/',sep='')
dir.create(pca_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pc <- function( df, path, file, var.axes ) {
    png(paste(path,file,sep=''), units="in", width=18.5, height=18.5, res=300)
    pca <- prcomp(t(df), center=TRUE, scale=TRUE)
    g <- ggbiplot(pca, obs.scale = 1, var.scale = 1,
    groups=f, ellipse=TRUE, circle=TRUE, var.axes=var.axes, varname.abbrev = FALSE)
    g <- g + scale_color_discrete(name = '')
    g <- g + theme(legend.direction = 'horizontal',legend.position = 'top')
    dev.off() }
pc(data, pca_path, "all_identified_pca.png", FALSE)
pc(data, pca_path, "all_identified_pca_labelled.png", TRUE)

if (length(rownames(pval_data)) >=2) {
    pc(pval_data, pca_path, "qval_0_05_identified_pca_labelled.png", TRUE)
    pc(pval_data, pca_path, "qval_0_05_identified_pca.png", FALSE)
}


# PCA significant
#pca_path=paste(outdir,'PCA/conditions/',sep='')
#dir.create(pca_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
#pc(cnd_pval_data, pca_path, "cnd_p_0_05_pca.png", TRUE)
#pc(cnd_pval_data_0_005, pca_path, "cnd_p_0_005_pca.png", TRUE)
#pc(cnd_pval_dup_data, pca_path, "cnd_p_0_05_mult_pca.png", TRUE)





#pc(data, pca_path, "all_identified_pca.png", FALSE)

#quit()

