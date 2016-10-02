#!/usr/bin/env Rscript

library("MSnID")
library("MSnbase")
library("dplyr")


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

#combined$accession <- combined$taxon_name
#combined["taxon_name"][is.na(combined[,"taxon_name"])] <- "No result"

peptide_map <- data.frame(exprs(msnset))
peptide_map <- add_rownames(peptide_map, "Row")
write.table(peptide_map, paste('peptide_map_pept2lca.txt',sep=''),sep='\t', row.names=FALSE)

psms(msnid) <- combined

msnset.count <- combineFeatures(msnset, fData(msnset)$taxon_name, redundancy.handler='unique', fun="sum",cv=FALSE)

peptide_counts <- data.frame(exprs(msnset))
peptide_counts <- add_rownames(peptide_counts, "Row")
write.table(peptide_counts, paste('peptide_sc_pept2lca.txt',sep=''),sep='\t', row.names=FALSE)

counts.df <- data.frame(exprs(msnset.count))
counts.df <- add_rownames(counts.df, "Row")
write.table(counts.df, paste('accession_sc_pept2lca.txt',sep=''),sep='\t', row.names=FALSE)

unipept.msnid <- msnid
save(unipept.msnid, file='msnid_processed_unipept_excl_decoy.Rdata')
msnset.taxon.count <- msnset.count
save(msnset.taxon.count, file='taxon_names_msnset.Rdata')


