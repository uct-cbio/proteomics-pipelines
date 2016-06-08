#!/opt/exp_soft/R-3.2.0/bin/Rscript
library("optparse")
library("BSgenome")

option_list = list(
  make_option(c("-g", "--genome"), type="character", default=NULL, 
              help="Genome fasta file, package will be generated in the same directory", metavar="character"), 
  make_option(c("-p", "--provider"), type="character", default=NULL, 
              help="Provider of the genome ie. NCBI, ENSEMBL", metavar="character"),
  make_option(c("-s", "--source_url"), type="character", default=NULL, 
              help="The permanent URL where the sequence data files used to forge the target package can be found.", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="The name of the genome assembly.", metavar="character"),
  make_option(c("-v", "--version"), type="character", default=NULL, 
              help="The version of the assembly.", metavar="character")

); 

opt_parser = OptionParser(option_list=option_list);
opt=parse_args(opt_parser);

genome <- opt$genome
provider <- opt$provider
source_url <- opt$source_url
name <- opt$name
version <- opt$version


con <- file(genome, open='r')
line <- strsplit(readLines(con)[[1]], split=',')[[1]][1]
words <- strsplit(line, split=' ')[[1]]

accession <- words[1]; words <- words[-1]
accession <- strsplit(accession, split='>')[[1]][2]

genus <- words[1]; words <- words[-1]
species <- words[1]; words <- words[-1]
strain <- paste(words,collapse='_')

abbr_org <- paste(genus, species, strain, sep ='_')
#abbr_org <- accession

split <- strsplit(genome, split='/')[[1]]
genome_file <- tail(split, n=1)
length(split) <- length(split) -1
genome_path = paste(split,collapse='/')
print(genome_file)
print(genome_path)

package <- paste('BSgenome', name, provider, version, sep='.')
title <- paste('Full genome sequence for', genus, species, strain,'-', accession, sep=' ')
description <- line
version <- '1.0.0'
author <- provider
maintainer <- 'name surname <email>'
license <- provider
suggests <- 'None'
organism <- paste(genus, species, strain, sep =' ')
common_name <- organism
provider <- provider
provider_version <- accession
release_date <- 'None'
release_name <- accession
source_url <- source_url
organism_biocview <- paste(genus, species, strain, sep='_')
BSgenomeObjname <- abbr_org
seqs_srcdir <- genome_path
seqfile_name <- genome_file


seed_path <- paste(genome_path,"/", name, ".seed", sep='')

sink(seed_path)

cat(sprintf(paste('Package: ',package,sep='')), "\n")
cat(sprintf(paste('Title: ', title, sep='')), "\n")
cat(sprintf(paste('Description: ',description)),"\n")
cat(sprintf(paste('Version: ',version)),"\n")
cat(sprintf(paste('Author: ', author)),"\n")
cat(sprintf(paste('Maintainer: ', maintainer)),"\n")
cat(sprintf(paste('License: ',license)),'\n')
cat(sprintf(paste('Suggests: ',suggests)),'\n')
cat(sprintf(paste('organism: ',organism)),'\n')
cat(sprintf(paste('common_name: ',common_name)),'\n')
cat(sprintf(paste('provider: ',provider)),'\n')
cat(sprintf(paste('provider_version: ',provider_version)),'\n')
cat(sprintf(paste('release_date: ', release_date)),'\n')
cat(sprintf(paste('release_name: ',release_name)),'\n')
cat(sprintf(paste('source_url: ', source_url)),'\n')
cat(sprintf(paste('organism_biocview: ',organism_biocview)),'\n')
cat(sprintf(paste('BSgenomeObjname: ',BSgenomeObjname)),'\n')
cat(sprintf('seqnames: paste("chr",1, sep="")'),'\n')
cat(sprintf(paste('seqs_srcdir: ',seqs_srcdir)),'\n')

sink()
forgeBSgenomeDataPkg(seed_path)


