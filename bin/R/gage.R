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
make_option(c("-i", "--indir"),
       type="character",
       default=NULL,
       help="Directory containing geneset annotations",
        metavar="character"),
make_option(c("-k", "--keggid"),
        type="character",
        default=NULL,
        help="Directory containing MQ proteogenomics pipeline output files",
        metavar="character"),
make_option(c("-d", "--design"),
        type="character",
        default=NULL,
        help="Experimental design",
        metavar="character"),
make_option(c("-t", "--table"),
        type="character",
        default=NULL,
        help="Table with intensities",
        metavar="character"),
make_option(c("-g", "--genecol"),
        type="character",
        default=NULL,
        help="Table with intensities",
        metavar="character"),
make_option(c("-ko", "--kocol"),
        type="character",
        default="Leading.Protein.Kegg.Orthology.ID",
        help="Table with intensities",
        metavar="character"),
make_option(c("-p", "--pval"),
        type="double",
        default=0.05,
        help="P value cutoff",
        metavar="double")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

inpath = opt$indir
path = opt$outdir
species = opt$keggid
genecol = opt$genecol
pval = opt$pval
kocol =opt$kocol

source(opt$design)

# GENE SETS 
bp <- 'false'
bp_path <- paste(inpath,'/bpset.Rdata',sep='')
if(file.exists(bp_path)){
    load(bp_path)
    bp <- 'true'
}

mf <- 'false'
mf_path <- paste(inpath,'/mfset.Rdata',sep='')
if(file.exists(mf_path)){
    load(mf_path)
    mf <- 'true'
}
cc <- 'false'
cc_path <- paste(inpath,'/ccset.Rdata',sep='')
if(file.exists(cc_path)){
    load(cc_path)
    cc <- 'true'
}
kegg <- 'false'
kegg_path <- paste(inpath,'/keggset.Rdata',sep='')
if(file.exists(kegg_path)){
    load(kegg_path)
    kegg <- 'true'
}
ec <- 'false'
ec_path <- paste(inpath,'/ecset.Rdata',sep='')
if(file.exists(ec_path)){
    load(ec_path)
    ec <- 'true'
}
operon <-'false'
operon_path <- paste(inpath,'/operonset.Rdata',sep='')
if(file.exists(operon_path)){
    load(operon_path)
    operon <- 'true'
}

ipr <- 'false'
ipr_path <- paste(inpath,'/iprset.Rdata',sep='')
if(file.exists(ipr_path)){
    load(ipr_path)
    ipr <- 'true'
}

metacyc <- 'false'
metacyc_path <- paste(inpath,'/metacycset.Rdata',sep='')
if(file.exists(metacyc_path)){
    load(metacyc_path)
    metacyc <- 'true'
}

reactome <- 'false'
reactome_path <- paste(inpath,'/reactomeset.Rdata',sep='')
if(file.exists(reactome_path)){
    load(reactome_path)
    reactome <- 'true'
}

table_path <- opt$table 
infile <- basename(table_path)

table <- read.csv(table_path)
row.names(table) <- table$Row.names

table[, cols] <- lapply(table[, cols], function(x){replace(x, x == 0,  NA)})
table[, cols] <- lapply(table[, cols], function(x){ log2(x)})


# GI table
gi_table <- cbind(table)
rownames(gi_table) <- gi_table$Identifier

#print(head(gi_table))
#quit()

s <- strsplit(as.character(gi_table[,genecol]), split = ";")
newtable <- data.frame( Identifier = rep(gi_table$Identifier, sapply(s, length)), gi = unlist(s))
refdata <- merge( table, newtable, by="Identifier", all.y = TRUE)
refdata <- refdata[!duplicated(refdata$gi), ]
ref <- refdata$gi
refdata <- refdata[,cols]
row.names(refdata) <-ref


ko_s <- strsplit(as.character(gi_table[,kocol]), split = ";")
ko_newtable <- data.frame( Identifier = rep(gi_table$Identifier, sapply(ko_s, length)), gi = unlist(ko_s))
kodata <- merge( table, ko_newtable, by="Identifier", all.y = TRUE)
kodata <- kodata[!duplicated(kodata$gi), ]
ko <- kodata$gi
kodata <- kodata[,cols]
row.names(kodata) <- ko


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

cutoff <- pval
print(cutoff)

less <- function(res, samp, ref) {
  less <- as.data.frame(res$less)
  #less <- less[!is.na(less$`p.val`),]
  less$Exposed <- samp
  less$Control <- ref
  less$RowName <- as.character(row.names(less))
  less$Coregulated <- "Down"
  if (length(row.names(less)) > 0) {
    less <- less[less$`p.val` < cutoff, ]
  }
  less <- less[!is.na(less$`p.val`),]
  #print(head(less))
  return(less)
}

greater <- function(res, samp, ref) {
  greater <- as.data.frame(res$greater)
  greater$RowName <- as.character(row.names(greater))
  greater$Exposed <- samp
  greater$Control <- ref
  greater$Coregulated <- "Up"
  if (length(row.names(greater)) > 0) {
    greater <- greater[greater$`p.val` < cutoff, ]
  }
  greater <- greater[!is.na(greater$`p.val`),]
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
    write.csv(gt, paste('IPR.up.', infile, sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls, paste('IPR.down.', infile, sep=''))
  }
 

  # EC
  print("EC")
  #operon_table <- table[row.names(table) %in% operon.set,]
  res <- analyse(table, ec.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "True"
    write.csv(gt, paste('EC.up.', infile,sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls,paste('EC.down.', infile, sep=''))
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
  ref.d <- refdata[, sampcols]-rowMeans(refdata[, refcols, drop=FALSE])
  ko.d <- kodata[, sampcols]-rowMeans(kodata[, refcols, drop=FALSE])
  
  res <- analyse(table, kegg.set, refcols, sampcols, T)
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    ls$RowName = paste(species, ls$RowName,sep = "")
    write.csv(ls, paste('KEGG.down.', infile, sep=''))
    less_ids <- row.names(ls) 
    try(pv.out.list <- sapply(less_ids, function(pid) pathview(gene.data = ref.d, kegg.native = T, out.suffix = infile, same.layer = F, pathway.id = paste(species, pid, sep=''), species = species)))
    try(pv.out.list <- sapply(less_ids, function(pid) pathview(gene.data = ko.d, kegg.native = T, out.suffix = infile, same.layer = F, pathway.id = paste('ko', pid, sep=''), species = 'ko')))
  }
  
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "True"
    gt$RowName = paste(species, gt$RowName,sep = "")
    write.csv(gt, paste('KEGG.up.', infile, sep=''))
    greater_ids <- row.names(gt)
    try(pv.out.list <- sapply(greater_ids, function(pid) pathview(gene.data = ref.d, kegg.native = T, out.suffix = infile, same.layer = F, pathway.id = paste(species, pid, sep=''), species = species)))
    try(pv.out.list <- sapply(greater_ids, function(pid) pathview(gene.data = ko.d, kegg.native = T, out.suffix = infile, same.layer = F, pathway.id = paste('ko', pid, sep=''), species = 'ko')))
  }
  
  res <- analyse(table, kegg.set, refcols, sampcols, F)
  gt <- greater(res, samp, ref)
  ls <- less(res, samp, ref)
  print(str(res))

  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "False"
    gt$RowName = paste(species, gt$RowName,sep = "")
    write.csv(gt, paste('KEGG.both.', infile,sep= ''))
    greater_ids <- row.names(gt)
  try(pv.out.list <- sapply(greater_ids, function(pid) pathview(gene.data = ref.d, kegg.native = T, out.suffix = infile, same.layer = F, pathway.id = paste(species, pid, sep=''), species = species)))
  try(pv.out.list <- sapply(greater_ids, function(pid) pathview(gene.data = ko.d, kegg.native = T, out.suffix = infile, same.layer = F, pathway.id = paste('ko', pid, sep=''), species = 'ko')))
  } 
  
  # BP
  print("BP")
  res <- analyse(table, bp.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
      gt$SameDir <- "True"
      write.csv(gt,paste('BP.up.',infile, sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls, paste('BP.down.', infile, sep=''))
  }
 
  # MF
  print("MF")
  res <- analyse(table, mf.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "True"
    write.csv(gt, paste( 'MF.up.', infile, sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls, paste( 'MF.down.', infile,sep= ''))
  }
 
  # CC
  print("CC")
  res <- analyse(table, cc.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "True"
    write.csv(gt, paste('CC.up.', infile, sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls, paste('CC.down.',infile,sep= ''))
  }

  # OPERONS
  if (operon == 'true') {
  print("OPERONS")
  #operon_table <- table[row.names(table) %in% operon.set,]
  res <- analyse(table, operon.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "True"
    write.csv(gt, paste('OPERON.up.',infile, sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls, paste('OPERON.down.', infile, sep=''))
   }}
  # RECTOME
  if (reactome == 'true') {
  print("REACTOME")
  res <- analyse(table, reactome.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "True"
    write.csv(gt, paste('REACTOME.up.', infile,sep= ''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls,paste( 'REACTOME.down.', infile , sep=''))
   }}
  # METACYC
  if (metacyc == 'true') {
  print("METACYC")
  res <- analyse(table, metacyc.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    gt$SameDir <- "True"
    write.csv(gt, paste( 'METACYC.up.', infile, sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    ls$SameDir <- "True"
    write.csv(ls, paste('METACYC.down.', infile, sep=''))
   }}
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

    outpath <-paste( path, '/', samp, '_', ref ,sep='' ) 
    process(table, refcols, sampcols, outpath, refdata, samp, ref) }
