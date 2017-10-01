#!/usr/bin/env Rscript

library('gage')
library('limma')
library('pathview')
library('optparse')


option_list = list(
make_option(c("-o", "--outdir"),
       type="character",
       default=NULL,
       help="Directory containing MQ proteogenomics pipeline output files",
        metavar="character"),
make_option(c("-k", "--keggid"),
        type="character",
        default=NULL,
        help="Directory containing MQ proteogenomics pipeline output files",
        metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

path = opt$outdir
species = opt$keggid
     
source(paste( path, '/experimental_design.R',sep=''))

# GENE SETS 
load(paste(path,'/gsea/bpset.Rdata',sep=''))
load(paste(path,'/gsea/mfset.Rdata',sep=''))
load(paste(path,'/gsea/ccset.Rdata',sep=''))
load(paste(path,'/gsea/keggset.Rdata',sep=''))
load(paste(path,'/gsea/operonset.Rdata',sep=''))
load(paste(path,'/gsea/iprset.Rdata',sep=''))

table_path <- paste(path, '/diff/msnbase/imputed.csv', sep='')
table <- read.csv(table_path)
row.names(table) <- table$Row.names



# GI table
gi_table <- cbind(table)

rownames(gi_table) <- gi_table$Identifier

#print(head(gi_table))
#quit()

s <- strsplit(as.character(gi_table$GeneID), split = ";")

newtable <- data.frame( Identifier = rep(gi_table$Identifier, sapply(s, length)), gi = unlist(s))

refdata <- merge( table, newtable, by="Identifier", all.y = TRUE)

#refdata$gi <- as.numeric(refdata$gi)

refdata <- refdata[!duplicated(refdata$gi), ]

ref <- refdata$gi
refdata <- refdata[,cols]

row.names(refdata) <-ref
#print(head((refdata))

#print(rownames(refdata))

ids <- table$Identifier

table <- table[,cols]

row.names(table) <- ids


#outpath <-paste(path,'/gsea/comparisons',sep='')
#dir.create(outpath, showWarnings = FALSE)

analyse <- function(data, gset, refcols, sampcols, samedir) {
  cnts.p <- gage(data, gsets = gset, ref = refcols, samp = sampcols, compare ="unpaired", samedir= samedir)  
  return(cnts.p)
}

less <- function(res, samp, ref) {
  less <- as.data.frame(res$less)
  less <- less[!is.na(less$`p.val`),]
  less$Exposed <- samp
  less$Control <- ref
  less$RowName <- as.character(row.names(less))
  less$Coregulated <- "Down"
  if (length(row.names(less)) > 0) {
    less <- less[less$`p.val` <= 0.05, ]
  }
  #print(head(less))
  return(less)
}

greater <- function(res, samp, ref) {
  greater <- as.data.frame(res$greater)
  greater <- greater[!is.na(greater$`p.val`),]
  greater$RowName <- as.character(row.names(greater))
  greater$Exposed <- samp
  greater$Control <- ref
  greater$Coregulated <- "Up"
  if (length(row.names(greater)) > 0) {
    greater <- greater[greater$`p.val` <= 0.05, ]
  }
  #print(head(greater))
  return(greater)
}

process <- function(table , refcols, sampcols, outpath, refdata, samp, ref) {
  # IPR
  print("IPR")
  
  dir.create(outpath, showWarnings = FALSE)
  setwd(outpath)
  #operon_table <- table[row.names(table) %in% operon.set,]
  res <- analyse(table, ipr.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "True"
    write.csv(gt, paste(outpath, '/IPR.up.csv', sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls, paste(outpath, '/IPR.down.csv',sep=''))
  }
  
  #res <- analyse(table, ipr.set, refcols, sampcols, FALSE)
  #gt <- greater(res, samp, ref)
  #if (length(row.names(gt)) > 0) {
  #  gt$SameDir <- "False"
  #  write.csv(gt, paste(outpath, '/IPR.up.both.csv', sep=''))
  #}
  #ls <- less(res, samp, ref)
  #if (length(row.names(ls)) > 0) {
  #  ls$SameDir <- "False"
  #  write.csv(ls, paste(outpath, '/IPR.down.both.csv',sep=''))
  #}




  print("KEGG")
  
  ref.d <- refdata[, sampcols]-rowMeans(refdata[, refcols])
  
  res <- analyse(table, kegg.set, refcols, sampcols, TRUE)

  ls <- less(res, samp, ref)

  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    
    ls$RowName = paste(species, ls$RowName,sep = "")

    write.csv(ls, paste(outpath, '/KEGG.down.csv',sep=''))
    less_ids <- row.names(ls) 
    pv.out.list <- sapply(less_ids, function(pid) pathview(gene.data = ref.d, kegg.native = T, out.suffix = 'keggViewDownregulated', same.layer = F, pathway.id = paste(species, pid, sep=''), species = species))
    pv.out.list <- sapply(less_ids, function(pid) pathview(gene.data = ref.d, kegg.native = F, out.suffix = 'graphvizViewDownregulated', split.group = T, same.layer = F, pathway.id = paste(species, pid, sep=''), species = species))
  }
  
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "True"
    gt$RowName = paste(species, gt$RowName,sep = "")
    write.csv(gt, paste(outpath, '/KEGG.up.csv',sep=''))
    greater_ids <- row.names(gt)
    pv.out.list <- sapply(greater_ids, function(pid) pathview(gene.data = ref.d, kegg.native = T, out.suffix = 'keggViewUpregulated', same.layer = F, pathway.id = paste(species, pid, sep=''), species = species))
    pv.out.list <- sapply(greater_ids, function(pid) pathview(gene.data = ref.d, kegg.native = F, out.suffix = 'graphvizViewUpregulated', split.group = T, same.layer = F, pathway.id = paste(species, pid, sep=''), species = species))
  }
  
  
  res <- analyse(table, kegg.set, refcols, sampcols, FALSE)
  gt <- greater(res, samp, ref)
  ls <- less(res, samp, ref)
  print(str(res))

  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "False"
    gt$RowName = paste(species, gt$RowName,sep = "")
    write.csv(gt, paste(outpath, '/KEGG.up.both.csv',sep=''))
    greater_ids <- row.names(gt)
  }
  #pv.out.list <- sapply(greater_ids, function(pid) pathview(gene.data = ref.d, kegg.native = T, out.suffix = 'keggViewUpBoth', same.layer = F, pathway.id = paste(species, pid, sep=''), species = species))
  #pv.out.list <- sapply(greater_ids, function(pid) pathview(gene.data = ref.d, kegg.native = F, out.suffix = 'graphvizViewUpBoth', split.group = T, same.layer = F, pathway.id = paste(species, pid, sep=''), species = species))
  
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "False"
    ls$RowName = paste(species, ls$RowName,sep = "")
    write.csv(ls, paste(outpath, '/KEGG.down.both.csv',sep=''))
    lesser_ids <- row.names(ls)
  }
  
  #pv.out.list <- sapply(lesser_ids, function(pid) pathview(gene.data = ref.d, kegg.native = T, out.suffix = 'keggViewDownBoth', same.layer = F, pathway.id = paste(species, pid, sep=''), species = species))
  #pv.out.list <- sapply(lesser_ids, function(pid) pathview(gene.data = ref.d, kegg.native = F, out.suffix = 'graphvizViewDownBoth', split.group = T, same.layer = F, pathway.id = paste(species, pid, sep=''), species = species))
  
  # BP
  print("BP")
  res <- analyse(table, bp.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
      gt$SameDir <- "True"
      write.csv(gt, paste(outpath, '/BP.up.csv',sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls, paste(outpath, '/BP.down.csv',sep=''))
  }
 
  # MF
  print("MF")
  res <- analyse(table, mf.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "True"
    write.csv(gt, paste(outpath, '/MF.up.csv',sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls, paste(outpath, '/MF.down.csv',sep=''))
  }
 
  # CC
  print("CC")
  res <- analyse(table, cc.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "True"
    write.csv(gt, paste(outpath, '/CC.up.csv', sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls, paste(outpath, '/CC.down.csv',sep=''))
  }

  # OPERONS
  print("OPERONS")
  #operon_table <- table[row.names(table) %in% operon.set,]
  res <- analyse(table, operon.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "True"
    write.csv(gt, paste(outpath, '/OPERON.up.csv', sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls, paste(outpath, '/OPERON.down.csv',sep=''))
   }
}

# Samples 
cols <- cols # Defined in experimental design template
f <- f # Defined in experimental design template
refmap <- data.frame(f, cols)
comparisons <- colnames(contrast.matrix)

for ( comp in comparisons){
    vals <- strsplit(comp,'-')
    # Get the reference cols
    ref  = as.character(vals[[1]][2])
    refdf <- refmap[refmap$f==ref,]
    #refcols <- as.character(refdf$cols)
    refcols <- as.numeric(rownames(refdf))
    
    # Get the sample cols
    samp = as.character(vals[[1]][1])
    sampdf <- refmap[refmap$f==samp,]
    sampcols <- as.numeric(rownames(sampdf))

    outpath <-paste( path, '/gsea/', samp, '_', ref ,sep='' ) 
    process(table, refcols, sampcols, outpath, refdata, samp, ref) }
