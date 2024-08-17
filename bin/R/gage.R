#!/usr/bin/env Rscript

library('gage')
library('limma')
library('pathview')
library('optparse')
library('tools')

print("Starting gene set enrichment")
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
make_option(c("-c", "--kocol"),
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
table_path <- opt$table 

source(opt$design)

# GENE SETS 
bp <- FALSE
bp_path <- paste(inpath,'/bpset.Rdata',sep='')
if(file.exists(bp_path)){
    load(bp_path)
    bp <- TRUE
}

mf <- FALSE
mf_path <- paste(inpath,'/mfset.Rdata',sep='')
if(file.exists(mf_path)){
    load(mf_path)
    mf <- TRUE
}
cc <- FALSE
cc_path <- paste(inpath,'/ccset.Rdata',sep='')
if(file.exists(cc_path)){
    load(cc_path)
    cc <- TRUE
}
kegg <- FALSE
kegg_path <- paste(inpath,'/keggset.Rdata',sep='')
if(file.exists(kegg_path)){
    load(kegg_path)

    kegg <- TRUE
}

kopathway <- FALSE
kopathway_path <- paste(inpath,'/kopathwayset.Rdata',sep='')
if(file.exists(kopathway_path)){
    load(kopathway_path)
    kopathway <- TRUE
}


ec <- FALSE
ec_path <- paste(inpath,'/ecset.Rdata',sep='')
if(file.exists(ec_path)){
    load(ec_path)
    ec <- TRUE
}
operon <-FALSE
operon_path <- paste(inpath,'/operonset.Rdata',sep='')
if(file.exists(operon_path)){
    load(operon_path)
    operon <- TRUE
}

ipr <- FALSE
ipr_path <- paste(inpath,'/iprset.Rdata',sep='')
if(file.exists(ipr_path)){
    load(ipr_path)
    ipr <- TRUE
}

metacyc <- FALSE
metacyc_path <- paste(inpath,'/metacycset.Rdata',sep='')
if(file.exists(metacyc_path)){
    load(metacyc_path)
    metacyc <- TRUE
}

reactome <- FALSE
reactome_path <- paste(inpath,'/reactomeset.Rdata',sep='')
if(file.exists(reactome_path)){
    load(reactome_path)
    reactome <- TRUE
}

infile <- basename(table_path)

table <- read.csv(table_path)

table[, cols] <- lapply(table[, cols], function(x){replace(x, x == 0,  NA)})
table[, cols] <- lapply(table[, cols], function(x){ log2(x)})


# GI table
gi_table <- cbind(table)
#ids <-  as.character(gi_table$Identifier)
#genecol <- as.character(gi_table[,genecol])
#gi_table <- gi_table[,cols]
#row.names(gi_table) <- ids
#gi_table
#quit()
s <- strsplit(as.character(gi_table[,genecol]), split = ";")
newtable <- data.frame( Identifier = rep(gi_table$Identifier, sapply(s, length)), gi = unlist(s))
refdata <- merge( table, newtable, by="Identifier", all.y = TRUE)

#
refdata$rowMean <- rowMeans(refdata[, cols])
refdata <- refdata[order(-refdata$rowMean),]


refdata <- refdata[!duplicated(refdata$gi), ]
ref <- refdata$gi
refdata <- refdata[,cols]
print('refcols')
#row.names(refdata) <-ref
.rowNamesDF(refdata, make.names=TRUE) <- ref

ko_s <- strsplit(as.character(gi_table[,kocol]), split = ";")
ko_newtable <- data.frame( Identifier = rep(gi_table$Identifier, sapply(ko_s, length)), gi = unlist(ko_s))
kodata <- merge( table, ko_newtable, by="Identifier", all.y = TRUE)

kodata$rowMean <- rowMeans(kodata[, cols])
kodata <- kodata[order(-kodata$rowMean),]
kodata <- kodata[!duplicated(kodata$gi), ]
ko <- kodata$gi
kodata <- kodata[,cols]
print('kocols')
.rowNamesDF(kodata, make.names=TRUE) <- ko




ids <- table$Identifier

table <- table[,cols]

print('idcols')
row.names(table) <- ids


#outpath <-paste(path,'/gsea/comparisons',sep='')
#dir.create(outpath, showWarnings = FALSE)

analyse <- function(data, gset, refcols, sampcols, samedir) {
  cnts.p <- gage(data, gsets = gset, ref = refcols, samp = sampcols, compare ="unpaired", same.dir= samedir) # used to be unpaired 
  return(cnts.p)
}

cutoff <- pval

current_wd = getwd()

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
  return(less)
}

get_species <- function(pid) {
    species <- gsub('[[:digit:]]+', '', pid)
    return(species)
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
  return(greater)
}

both <- function(res, samp, ref) {
  both <- as.data.frame(res$greater)
  #_ <- as.data.frame(res$less)
  both$RowName <- as.character(row.names(both))
  both$Exposed <- samp
  both$Control <- ref
  both$Coregulated <- "Both"
  if (length(row.names(both)) > 0) {
      both <- both[both$`p.val` < cutoff, ]
  }
  both <- both[!is.na(both$`p.val`),]
  return(both)
}

process <- function(table , refcols, sampcols, outpath, refdata, samp, ref) {
  dir.create(outpath, showWarnings = FALSE)
  custom_outpath <- paste(outpath, '/custom/', sep='')
  dir.create(custom_outpath, showWarnings = FALSE)
  setwd(outpath)
  #operon_table <- table[row.names(table) %in% operon.set,]
  if (ipr==TRUE) {
  # IPR
  print("IPR")
  res <- analyse(table, ipr.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    write.csv(gt, paste('IPR.up.', infile, sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    write.csv(ls, paste('IPR.down.', infile, sep=''))
  }}
  # EC
  #print("EC")
  ##operon_table <- table[row.names(table) %in% operon.set,]
  #res <- analyse(table, ec.set, refcols, sampcols, TRUE)
  #gt <- greater(res, samp, ref)
  #if (length(row.names(gt)) > 0) {
  #  gt$SameDir <- "True"
  #  write.csv(gt, paste('EC.up.', infile,sep=''))
  #}
  #ls <- less(res, samp, ref)
  #if (length(row.names(ls)) > 0) {
  #  ls$SameDir <- "True"
  #  write.csv(ls,paste('EC.down.', infile, sep=''))
  #}


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
  

  if (kegg==TRUE) {
  print("KEGG")
  ref.d <- as.matrix(refdata[, sampcols]-rowMeans(refdata[, refcols, drop=FALSE]))

  res <- analyse(table, kegg.set, refcols, sampcols, TRUE)
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    write.csv(ls, paste('KEGG.down.', infile, sep=''))
    less_ids <- row.names(ls) 
    try(pv.out.list <- sapply(less_ids, function(pid) pathview(gene.data = ref.d, kegg.native = T, out.suffix = infile, same.layer = F, pathway.id = pid, species = get_species(pid))))
  }
  
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    write.csv(gt, paste('KEGG.up.', infile, sep=''))
    greater_ids <- row.names(gt)
    try(pv.out.list <- sapply(greater_ids, function(pid) pathview(gene.data = ref.d, kegg.native = T, out.suffix = infile, gene.idtype='entrez', same.layer = F, pathway.id = pid, species = get_species(pid))))
  }
  
  res <- analyse(table, kegg.set, refcols, sampcols, FALSE)
  #print(summary(res))
  bt <- both(res, samp, ref)
  if (length(row.names(bt)) > 0) {
    write.csv(bt, paste('KEGG.both.', infile,sep= ''))
    both_ids <- row.names(bt)
  try(pv.out.list <- sapply(both_ids, function(pid) pathview(gene.data = ref.d, gene.idtype='entrez', kegg.native = T, out.suffix = infile, same.layer = F, pathway.id = pid, species = get_species(pid))))}}

  if (kopathway == TRUE) {
  print("KEGG KO PATHWAY")
  ko.d <- as.matrix(kodata[, sampcols] -rowMeans(kodata[, refcols, drop=FALSE]))
  res <- analyse(table, keggko.set, refcols, sampcols, TRUE)
  
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    write.csv(ls, paste('KOPATHWAY.down.', infile, sep=''))
    less_ids <- row.names(ls) 
    try(pv.out.list <- sapply(less_ids, function(pid) pathview(gene.data = ko.d, kegg.native = T, out.suffix = infile, same.layer = F, pathway.id = pid ,gene.idtype="kegg", species=get_species(pid))))
  }
  
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    #gt$RowName = paste(species, gt$RowName,sep = "")
    write.csv(gt, paste('KOPATHWAY.up.', infile, sep=''))
    greater_ids <- row.names(gt)
    try(pv.out.list <- sapply(greater_ids, function(pid) pathview(gene.data = ko.d, kegg.native = T, out.suffix = infile, same.layer = F, pathway.id = pid,gene.idtype="kegg", species=get_species(pid))))
  }
  
  res <- analyse(table, keggko.set, refcols, sampcols, FALSE)
  bt <- both(res, samp, ref)
  if (length(row.names(bt)) > 0) {
    write.csv(bt, paste('KOPATHWAY.both.', infile,sep= ''))
    both_ids <- row.names(bt)
    try(pv.out.list <- sapply(both_ids, function(pid) pathview(gene.data = ko.d, kegg.native = T, out.suffix = infile, same.layer = F, pathway.id = pid,gene.idtype="kegg", species=get_species(pid))))
  } }

  # BP
  if (bp == TRUE) {
  print("BP")
  
  res <- analyse(table, bp.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
      write.csv(gt,paste('BP.up.',infile, sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    write.csv(ls, paste('BP.down.', infile, sep=''))
  }}
 
  # MF
  if (mf == TRUE) {
  print("MF")
  res <- analyse(table, mf.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    write.csv(gt, paste( 'MF.up.', infile, sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    write.csv(ls, paste( 'MF.down.', infile,sep= ''))
  }}
 
  # CC
  if (cc == TRUE ) {
  print("CC")
  res <- analyse(table, cc.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    write.csv(gt, paste('CC.up.', infile, sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    write.csv(ls, paste('CC.down.',infile,sep= ''))
  }}

  # OPERONS
  if (operon == TRUE) {
  print("OPERONS")
  #operon_table <- table[row.names(table) %in% operon.set,]
  res <- analyse(table, operon.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    write.csv(gt, paste('OPERON.up.',infile, sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    write.csv(ls, paste('OPERON.down.', infile, sep=''))
   }}
  
  # RECTOME
  if (reactome == TRUE) {
  print("REACTOME")
  res <- analyse(table, reactome.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    write.csv(gt, paste('REACTOME.up.', infile,sep= ''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    write.csv(ls,paste( 'REACTOME.down.', infile , sep=''))
   }}
  # METACYC
  if (metacyc == TRUE) {
  print("METACYC")
  res <- analyse(table, metacyc.set, refcols, sampcols, TRUE)
  gt <- greater(res, samp, ref)
  if (length(row.names(gt)) > 0) {
    write.csv(gt, paste( 'METACYC.up.', infile, sep=''))
  }
  ls <- less(res, samp, ref)
  if (length(row.names(ls)) > 0) {
    write.csv(ls, paste('METACYC.down.', infile, sep=''))
   }}
   
  # CUSTOM GENESETS
  setwd(current_wd)
  custom_inpath <- paste(inpath,'/custom/', sep='')
  files <- list.files(custom_inpath)

  for (file in files)
  {
      if (file_ext(file) == 'Rdata') {
          
          
          load( paste(custom_inpath,file, sep=''))
          
          setwd(custom_outpath)
          filename <- sub('\\.csv.Rdata$', '', file) 
          print(paste("PROCESSING: ", filename, sep=''))
          
          res <- analyse(table, gtab.set, refcols, sampcols, TRUE)
          gt <- greater(res, samp, ref)
          if (length(row.names(gt)) > 0) {
          prefix <- paste(filename,'.up.',sep='')
          write.csv(gt, paste( prefix, infile, sep=''))
          }
          ls <- less(res, samp, ref)
          if (length(row.names(ls)) > 0) {
            prefix <- paste(filename,'.down.',sep='')
            write.csv(ls, paste(prefix, infile, sep=''))
          }
          setwd(current_wd)

      }
  }

}

# Samples 
cols <- cols # Defined in experimental design template
f <- f # Defined in experimental design template
refmap <- data.frame(f, cols)
comparisons <- colnames(contrast.matrix)



for ( comp in comparisons){
    setwd(current_wd)
    #print("here")
    #print(comp)
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
    #print(head(table))
    #print(refcols)
    #print(sampcols)

    process(table, refcols, sampcols, outpath, refdata, samp, ref) }
