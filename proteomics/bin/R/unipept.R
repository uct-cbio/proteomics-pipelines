#!/usr/bin/env Rscript

library("MSnID")
library("MSnbase")
library("dplyr")

options(bitmapType='cairo') 

msnid <- MSnID(".")

p2lca <- read.csv('pept2lca.csv')
load('msnid_processed.Rdata')
processed.msnid <- apply_filter(processed.msnid, "isDecoy == FALSE")
ids <- psms(processed.msnid)
msnset <- as(processed.msnid, "MSnSet")
fdata <- fData(msnset)

fdata$pepSeq <- ids$peptideSequence[match(fdata$peptide, ids$peptide)]
fdata$taxon_name <- p2lca$taxon_name[match(fdata$pepSeq, p2lca$peptide)]

fdata$taxon_name <- as.character(fdata$taxon_name)
fdata$taxon_name[is.na(fdata$taxon_name)] <- "unassigned"


fdata$taxon_name=as.factor(fdata$taxon_name)
fData(msnset) <- fdata
combined <- merge(ids, p2lca, by.x='peptideSequence', by.y='peptide', all.x=TRUE)
write.table(combined, paste('merged_msnid_unipept.txt',sep=''),sep='\t', row.names=FALSE)

psms(msnid) <- combined

msnset.count <- combineFeatures(msnset, fData(msnset)$taxon_name, redundancy.handler='unique', fun="sum",cv=FALSE)

peptide_counts <- data.frame(exprs(msnset))
peptide_counts <- add_rownames(peptide_counts, "Row")
write.table(peptide_counts, paste('peptide_sc_pept2lca.txt',sep=''),sep='\t', row.names=FALSE)

counts.df <- data.frame(exprs(msnset.count))
counts.df$TotalMSMS <- rowSums(counts.df)
counts.df <- add_rownames(counts.df, "Row")
counts.df <- counts.df[with(counts.df, order(-TotalMSMS)), ]

write.table(counts.df, paste('accession_sc_pept2lca.txt',sep=''),sep='\t', row.names=FALSE)
# Create PIE
# PSM condifence


jpeg('accession_sc_pie.jpeg', width=1000,height=900)
filt <- counts.df[!(counts.df$Row=="unassigned"),]
#par(mar=c(6,12,6,12)+.1)
slices <- round((filt$TotalMSMS/sum(filt$TotalMSMS) * 100),3)
lbls <- paste(filt$Row, " ", slices,' %', sep="")
df <- data.frame(matrix(, nrow=length(lbls), ncol=2) )
df$label <- lbls
df$slices <- slices
df$label[df$slices < 0.5] <- ""
par(mar=c(6,12,6,12)+.1)
pie(df$slices, labels =df$label ,  main="Pie Chart of total spectral counts by Unipept LCA (%)")
dev.off()

unipept.msnid <- msnid
save(unipept.msnid, file='msnid_processed_unipept_excl_decoy.Rdata')
msnset.taxon.count <- msnset.count
save(msnset.taxon.count, file='taxon_names_msnset.Rdata')


