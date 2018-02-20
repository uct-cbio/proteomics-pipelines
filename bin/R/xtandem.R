#!/usr/bin/env Rscript
library("optparse")
library(rTANDEM)

options(bitmapType='cairo')

option_list = list(
make_option(c("-m", "--mgf"), type="character", default=NULL,
              help="Directory of mzIdentML files for global FDR control", metavar="character"),
make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="Percentage value to control the FDR", metavar="character"),
make_option(c("-ixml","--input_xml"), type="character", default='NULL',
              help="default_input.xml", metavar="character"),
make_option(c("-ixsl","--input_xsl"), type="character", default='NULL',
              help="tandem-input-style.xsl", metavar="character"),
make_option(c("-oxml","--output_xml"), type="character", default='NULL',
              help="output.xml", metavar="character"))

opt_parser = OptionParser(option_list=option_list);

opt = parse_args(opt_parser);

mgf = opt$mgf
fasta = opt$fasta
input_xml = opt$input_xml
input_xsl = opt$input_xsl
output_xml = opt$output_xml


print(mgf)
print(fasta)

taxonomy <- rTTaxo( taxon="metanovo", format="peptide", fasta)
param <- rTParam()
param <- setParamValue(param, 'protein', 'taxon', value="metanovo")
param <- setParamValue(param, 'list path', 'taxonomy information', taxonomy)
param <- setParamValue(param, 'list path', 'default parameters', input_xml)

param <- setParamValue(param, 'spectrum', 'path', mgf)
param <- setParamValue(param, 'output', 'xsl path', input_xsl)
param <- setParamValue(param, 'output', 'path', output_xml)
result.path <- tandem(param)

