#!/usr/bin/env Rscript
#options(error=traceback)
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
    go <- GOTERM[goid]
    return(Term(go)[1] ) }

go.ontology <- function(goid ) {
    go <- GOTERM[goid]
    ontology <- Ontology(go)[1]
    return(ontology)}
gopath <- paste(outpath,'/go2proteingroups.csv', sep='') 
gtab <- read.csv(gopath)
gtab <- data.frame(lapply(gtab, as.character), stringsAsFactors=FALSE)
gtab$GO_TERM <- lapply(gtab$GO_ID, go.term)
gtab$GO_ONTOLOGY <- lapply(gtab$GO_ID, go.ontology)
gtab$ID<- with(gtab, paste0(GO_ID, " ", GO_TERM))
gtab$GENES <- lapply(gtab$GENES, vct)
# BIOLOGICAL PROCESS
bp <- gtab[gtab$GO_ONTOLOGY == 'BP',]
bp.set <- as.list(bp$GENES)
names(bp.set) <- as.list(bp$ID)
save(bp.set, file=paste(outpath, '/bpset.Rdata', sep=''))

# CELLULAR COMPONENT
cc <- gtab[gtab$GO_ONTOLOGY == 'CC',]
cc.set <- as.list(cc$GENES)
names(cc.set) <- as.list(cc$ID)
save(cc.set, file=paste(outpath, '/ccset.Rdata', sep=''))

# MOLECULAR FUNCTION
mf <- gtab[gtab$GO_ONTOLOGY == 'MF',]
mf.set <- as.list(mf$GENES)
names(mf.set) <- as.list(mf$ID)
save(mf.set, file=paste(outpath, '/mfset.Rdata', sep=''))

###############
#    KEGG     #
###############

keggpath <- paste(outpath,'/kegg2proteingroups.csv', sep='') 
ktab <- read.csv(keggpath, colClasses=c("character","character"))
ktab <- data.frame(lapply(ktab, as.character), stringsAsFactors=FALSE)

#ktab$ID<- with(ktab, paste0(kegg_id, KEGG_ID))

ktab$ID <- ktab$KEGG_ID
ktab$GENES <- lapply(ktab$GENES, vct)

kegg.set <- as.list(ktab$GENES)
names(kegg.set) <- as.list(ktab$ID)
save(kegg.set, file=paste(outpath, '/keggset.Rdata', sep=''))

###############
#    EC     #
###############

ecpath <- paste(outpath,'/ec2proteingroups.csv', sep='') 
etab <- read.csv(ecpath, colClasses=c("character","character"))
etab <- data.frame(lapply(etab, as.character), stringsAsFactors=FALSE)

#ktab$ID<- with(ktab, paste0(kegg_id, KEGG_ID))

etab$ID <- etab$EC_ID
etab$GENES <- lapply(etab$GENES, vct)

ec.set <- as.list(etab$GENES)
names(ec.set) <- as.list(etab$ID)
save(ec.set, file=paste(outpath, '/ecset.Rdata', sep=''))

###############
#    METACYC     #
###############

metacycpath <- paste(outpath,'/metacyc2proteingroups.csv', sep='') 
mtab <- read.csv(metacycpath, colClasses=c("character","character"))
mtab <- data.frame(lapply(mtab, as.character), stringsAsFactors=FALSE)
mtab$ID <- mtab$METACYC_ID
mtab$GENES <- lapply(mtab$GENES, vct)
metacyc.set <- as.list(mtab$GENES)
names(metacyc.set) <- as.list(mtab$ID)
save(metacyc.set, file=paste(outpath, '/metacycset.Rdata', sep=''))

###############
#    REACTOME    #
###############

reactomepath <- paste(outpath,'/reactome2proteingroups.csv', sep='') 
rtab <- read.csv(reactomepath, colClasses=c("character","character"))
rtab <- data.frame(lapply(rtab, as.character), stringsAsFactors=FALSE)
rtab$ID <- rtab$REACTOME_ID
rtab$GENES <- lapply(rtab$GENES, vct)
reactome.set <- as.list(rtab$GENES)
names(reactome.set) <- as.list(rtab$ID)
save(reactome.set, file=paste(outpath, '/reactomeset.Rdata', sep=''))

#################
#    OPERONS    #
#################

operonpath <- paste(outpath,'/operons2proteingroups.csv', sep='') 
if(file.exists(operonpath)){
    operontab <- read.csv(operonpath, colClasses=c("character","character"))
    operontab <- data.frame(lapply(operontab, as.character), stringsAsFactors=FALSE)

    #ktab$ID<- with(ktab, paste0(kegg_id, KEGG_ID))

    operontab$ID <- operontab$OPERON_ID
    operontab$GENES <- lapply(operontab$GENES, vct)

    operon.set <- as.list(operontab$GENES)
    names(operon.set) <- as.list(operontab$ID)
    save(operon.set, file=paste(outpath, '/operonset.Rdata', sep=''))

}



#################
#    IPR        #
#################
iprpath <- paste(outpath,'/ipr2proteingroups.csv', sep='') 
iprtab <- read.csv(iprpath, colClasses=c("character","character"))
iprtab <- data.frame(lapply(iprtab, as.character), stringsAsFactors=FALSE)
#ktab$ID<- with(ktab, paste0(kegg_id, KEGG_ID))
iprtab$ID <- iprtab$IPR_ID
iprtab$GENES <- lapply(iprtab$GENES, vct)
ipr.set <- as.list(iprtab$GENES)
names(ipr.set) <- as.list(iprtab$ID)
save(ipr.set, file=paste(outpath, '/iprset.Rdata', sep=''))
