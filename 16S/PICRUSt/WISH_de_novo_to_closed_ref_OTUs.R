#PURPOSE OF THIS SCRIPT: TO USE EXISTING DE NOVO OTU TABLE (WISH FEMALE) AND USING THE RESULTS FROM MAPPING DE NOVO PICKED OTUS (USING THEIR SECQUENCES)
#AND USEARCH-GLOBAL TO GREENGENES IDS. BELOW: IMPORT THE DE NOVO OTU TABLE INTO R AND REPLACE THE ONES WITH GREENGENES MATCHES (TYPICALLY THOSE WITH THE BEST
#GREENGENES ANNOTATIONS
#THE RESULTING OTU TABLE CAN NOW BE USED TO CONDUCT PICRUST ANALYSIS
library(phyloseq)
setwd("/Users/jameslennard/Documents/Academic/Postdoc/projects/WISH_female") #set working directory as appropriate
inDir <- paste0(getwd(),"/process/PICRUSt/v3_May2016_redo_new_OTU_table") #set input directory as appropriate
outDir <- inDir #set output directory as appropriate

#import mapping list of de novo IDs to GG IDs
o <- read.table(paste0(inDir,"/de_novo_repset_to_GG_13_8_map_forR.txt"), sep = "\t", header =T) #import .uc file (converted to .txt file) obtained during standard cbio 16S pipeline
head(o)#note: and * in the 2nd column mean there was no GG ID match
length(which(o[,2]=="*"))#number of OTUs with no GG matches
o.closed <- o[o[,2]!="*",] #subset to exclude OTUs with no GG match
dim(o.closed)#number of OTUs that remain
head(o.closed)
rownames(o.closed) <- o.closed[,1]
#import OTU table
load(paste0(getwd(),"/process/RData/v7_primer_stripped_WISH_V1.RData"))#Change as appropriate - previously prepared OTU table from R, import as .RData object 

#---------------------------------------------------
#FILTER 'JUNK' OTUS AND STANDARDIZE OTU TABLE
M.prune <- prune_samples(sample_sums(M.x) > 5000, M.x)#In this case the loaded OTU .RData object was names 'M.x'
total = median(sample_sums(M.prune))
standf = function(x, t=total) round(t * (x / sum(x)))
M.std = transform_sample_counts(M.prune, standf)
#mild filter, just to get rid of the riff-raff
M.f = filter_taxa(M.std,function(x) sum(x > 10) > (0.02*length(x)) | sum(x) > 0.001*total, TRUE)#select only OTUs where at least 2% of samples have more than 10 counts of that OTU
#OR the rowsum for that OTU > 0.1% of the total median count (for cases where the minority of samples may have high counts of a rare OTU)
ntaxa(M.f)#Number of taxa that pass filter
#write OTU IDs to file to filter .biom format OTU table for use in PICRUSt
keep.IDs = rownames(otu_table(M.f))
length(keep.IDs)#
write.table(keep.IDs,file = paste0(getwd(),"/process/PICRUSt/v3_May2016_redo_new_OTU_table/OTU_IDs_after_mild_filter.txt"), quote = FALSE, row.names = FALSE)
#change output file specification as appropriate
o.tab <- otu_table(M.f)
dim(o.tab)#
length(which(rownames(o.tab)%in%o[,1]))#all match?
#---------------------------------------------------
#filter otu table to exclude de novo otus with no GG OTU matches
o.tab <- o.tab[rownames(o.tab)%in%o.closed[,1],]
dim(o.tab)#
#now substitute de novo IDs with GG IDs
o.closed <- o.closed[rownames(o.tab),]
rownames(o.tab) <- o.closed[,2]
head(o.tab)
dim(o.tab)#
#problem: how to deal with duplicate GG IDs - e.g. more than one of my de novo IDs map to the same GG ID
length(unique(rownames(o.tab)))#
head(duplicated(rownames(o.tab)))
#------
#example
rownames(o.tab)[85]
# rownames(o.tab)[85]
# [1] "198893"
# > cat("Synch1464614218233544000\n");
temp <- which(rownames(o.tab)=="198893")#84 85
o.tab[temp,]
x = o.closed[o.closed[,2]=="128382",]
x
# x
#         De_novo_ID GG_13_8_ID
# OTU_434    OTU_434     128382
# OTU_44      OTU_44     128382
# > cat("Synch1464614273695909000\n");
#end example
tax_table(M.x)[rownames(x),]#Dialister
#------
#Will need to collapse OTUs with identical GG IDs somehow by summing the counts for those OTUs and replacing it with that total
#step 1. sort table by ID
o.t <- o.tab[ order(row.names(o.tab)), ]
head(rownames(o.t))#example output
# head(rownames(o.t))
# [1] "1000986" "1007750" "1028501" "1033552" "1035532" "1038074"
# > cat("Synch1464614289894893000\n");
#---
#step 2:
library(plyr)
t = data.frame(o.t)
t$ID <- rownames(o.t)
test=ddply(t, .(ID), function(x) colSums(x[,-dim(x)[2]], na.rm = TRUE))#break into chunks according to the column 'ID' and sum each chunk
head(test)
rownames(test) <- test$ID
dim(test)#
test <- test[,-1]
head(test)#tada!!!
#export new otu table with closed ref (GG 13.8) IDs
#NB NB NB - you have to add in "OTUId" in the IDs column otherwise you'll get an error when you're trying to convert to .biom
write.table(test, paste0(outDir,"/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline.txt"), sep = "\t", quote = FALSE,col.names=NA)



