#!/opt/exp_soft/R-3.2.0/bin/Rscript 

library("dplyr")
library("optparse")

option_list = list(
 make_option(c("-t", "--taxid"), type="integer", default=NULL, 
              help="The taxid of interest", metavar="character"),
	make_option(c("-o","--out_file"), type="character", default='accession',
              help="The name of the output file (full path)", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(UniProt.ws)

up <- UniProt.ws(taxId=opt$t)

keys <- keys(up, keytype="UNIPROTKB")
columns <- columns(up)

df <- select(up, keys=keys, columns=columns, keytype="UNIPROTKB")

df <- df[, colSums(is.na(df)) != nrow(df)]

write.csv(df, file=opt$o, sep=',')

