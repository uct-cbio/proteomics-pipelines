#!/usr/bin/env Rscript

library('optparse')
library('GO.db')
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
outpath = opt$outdir
kegg_id = opt$keggid

vct <- function(chars) { return(as.vector(strsplit(chars, '|', fixed=TRUE))[[1]]) }


###################
#  GO Genesets    #
###################
go.term <- function(goid)
{
    go <- GOTERM[[goid]]
    return(go@Term)  }

go.ontology <- function(goid)
{
    go <- GOTERM[[goid]]
    return(go@Ontology) }


gopath <- paste(outpath,'/gsea/go2proteingroups.csv', sep='') 
gtab <- read.csv(gopath)
gtab <- data.frame(lapply(gtab, as.character), stringsAsFactors=FALSE)
#gtab <- head(gtab)

gtab$GO_TERM <- lapply(gtab$GO_ID, go.term)
gtab$GO_ONTOLOGY <- lapply(gtab$GO_ID, go.ontology)
gtab$ID<- with(gtab, paste0(GO_ID, " ", GO_TERM))
gtab$GENES <- lapply(gtab$GENES, vct)


# BIOLOGICAL PROCESS
bp <- gtab[gtab$GO_ONTOLOGY == 'BP',]
bp.set <- as.list(bp$GENES)
names(bp.set) <- as.list(bp$ID)
save(bp.set, file=paste(outpath, '/gsea/bpset.Rdata', sep=''))

# CELLULAR COMPONENT
cc <- gtab[gtab$GO_ONTOLOGY == 'CC',]
cc.set <- as.list(cc$GENES)
names(cc.set) <- as.list(cc$ID)
save(cc.set, file=paste(outpath, '/gsea/ccset.Rdata', sep=''))

# MOLECULAR FUNCTION
mf <- gtab[gtab$GO_ONTOLOGY == 'MF',]
mf.set <- as.list(mf$GENES)
names(mf.set) <- as.list(mf$ID)
save(mf.set, file=paste(outpath, '/gsea/mfset.Rdata', sep=''))

###############
#    KEGG     #
###############

keggpath <- paste(outpath,'/gsea/kegg2proteingroups.csv', sep='') 
ktab <- read.csv(keggpath, colClasses=c("character","character"))
ktab <- data.frame(lapply(ktab, as.character), stringsAsFactors=FALSE)
#ktab$ID<- with(ktab, paste0(kegg_id, KEGG_ID))
ktab$ID <- ktab$KEGG_ID
ktab$GENES <- lapply(ktab$GENES, vct)

kegg.set <- as.list(ktab$GENES)
names(kegg.set) <- as.list(ktab$ID)
save(kegg.set, file=paste(outpath, '/gsea/keggset.Rdata', sep=''))

