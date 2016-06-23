#!/opt/exp_soft/R-3.2.0/bin/Rscript 

library("dplyr")
library("optparse")
option_list = list(
 make_option(c("-i", "--indir"), type="character", default=NULL, 
              help="Directory of mzIdentML files for gloabal FDR control", metavar="character"),
	make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Directory that will be created/overwritten with output files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$indir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input directory).n", call.=FALSE)
} 

unlink(paste(opt$indir,"/analysis",sep=''), recursive=TRUE)

dirfiles <- list.files(opt$indir, full.names=TRUE)
if (is.null(opt$outdir)) {
	outdir <- paste(opt$indir,'/analysis',sep='');
} else { outidr <- opt$outdir} 

dir.create(outdir,showWarnings=TRUE,recursive=FALSE,mode='0777')
sink(paste(outdir, '/log.txt',sep=''))

print(dirfiles)
library("MSnID")

msnid <- MSnID(".")
msnid <- read_mzIDs(msnid, dirfiles)
msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")

#all_descriptions <- msnid$description

jpeg(paste(outdir,'/missed_cleavages.jpeg',sep=''))
pepCleav <- unique(psms(msnid)[,c("numMissCleavages", "isDecoy", "peptide")])
pepCleav <- as.data.frame(table(pepCleav[,c("numMissCleavages", "isDecoy")]))
library("ggplot2")
ggplot(pepCleav, aes(x=numMissCleavages, y=Freq, fill=isDecoy)) + geom_bar(stat='identity', position='dodge')  + ggtitle("Number of Missed Cleavages")
dev.off()


jpeg(paste(outdir,'/peptide_lengths.jpeg',sep=''))
msnid$PepLength <- nchar(msnid$peptide) - 4
pepLen <- unique(psms(msnid)[,c("PepLength", "isDecoy", "peptide")])
ggplot(pepLen, aes(x=PepLength, fill=isDecoy)) + geom_histogram(position='dodge', binwidth=3) + ggtitle("Distribution of Peptide Lengths")
dev.off()

#parent ion mass measurement error (ppm)
jpeg(paste(outdir,'/parent_ion_mass_measurement_error.jpeg',sep=''))
ppm <- mass_measurement_error(msnid)
ggplot(as.data.frame(ppm), aes(x=ppm))  + geom_histogram(binwidth=100)
dev.off()


jpeg(paste(outdir,'/parent_ion_mass_measurement_targetdecoy.jpeg',sep=''))
dM <- with(psms(msnid),+ (experimentalMassToCharge-calculatedMassToCharge)*chargeState) 
x <- data.frame(dM, isDecoy=msnid$isDecoy)
ggplot(x, aes(x=dM, fill=isDecoy)) +geom_histogram(position='stack', binwidth=0.1)
dev.off()

# MSnID package provide a simple correct_peak_selection function that simply adds or subtracts the difference between 13C and 12C to make the error less then 1 Dalton.

jpeg(paste(outdir,'/parent_ion_mass_measurement_corrected.jpeg',sep=''))
msnid.fixed <- correct_peak_selection(msnid)
ppm <- mass_measurement_error(msnid.fixed)
ggplot(as.data.frame(ppm), aes(x=ppm)) + geom_histogram(binwidth=0.25)
dev.off()

# recalibrate
jpeg(paste(outdir,'/parent_ion_mass_measurement_corrected_recalibrated.jpeg',sep=''))
msnid <- recalibrate(msnid.fixed)
ppm <- mass_measurement_error(msnid)
ggplot(as.data.frame(ppm), aes(x=ppm)) +geom_histogram(binwidth=0.25)
dev.off()

# PSM condifence (converted to PEP scores) 
msnid$PeptideShakerPSMConfidence <- as.numeric(msnid$`PeptideShaker PSM confidence`)
msnid$PEP <- 1.0 - msnid$PeptideShakerPSMConfidence/100.0 

jpeg(paste(outdir,'/psm_PEP.jpeg',sep=''))
params <- psms(msnid)[,c("PEP","isDecoy")]
ggplot(params) + geom_density(aes(x = PEP, color = isDecoy, ..count..))
dev.off()

# PSM mass error

jpeg(paste(outdir,'/psm_absParentMassErrorPPM.jpeg',sep=''))
msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))
params <- psms(msnid)[,c("absParentMassErrorPPM","isDecoy")]
ggplot(params) + geom_density(aes(x = absParentMassErrorPPM, color = isDecoy, ..count..))
dev.off()

print('Creating MSnIDFilter object:')
filtObj <- MSnIDFilter(msnid)
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
filtObj$PEP <- list(comparison="<", threshold=0.01)
filtObj$PepLength <- list(comparison="<", threshold=30)
filtObj$numMissCleavages <- list(comparison="<", threshold=3)
filtObj$numIrregCleavages <- list(comparison="<", threshold=1)
show(filtObj)
evaluate_filter(msnid, filtObj, level='PSM')
cat('\n')

print('Optimizing MSnIDFilter object with "Grid" method:')
filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01,method="Grid",level="accession", n.iter=1000) 
show(filtObj.grid)
evaluate_filter(msnid, filtObj.grid, level='PSM')
cat('\n')

print('Optimizing MSnIDFilter.Grid object with "Nelder-Mead" method:')
filtObj.nm <- optimize_filter(filtObj.grid, msnid, fdr.max=0.01,method="Nelder-Mead",level="accession", n.iter=1000) 
show(filtObj.nm)
evaluate_filter(msnid, filtObj.nm, level='PSM')
cat('\n')

print('Optimizing MSnIDFilter.NM object with "Simulated-Annealing" method:')
filtObj.sa <- optimize_filter(filtObj.nm, msnid, fdr.max=0.01,method="SANN",level="accession", n.iter=1000) 
show(filtObj.sa)
evaluate_filter(msnid, filtObj.sa, level='PSM')
cat('\n')

#print(head(psms(msnid)))
#names(msnid)
#show(filtObj.sa)

print('Filter msnid object using the MSnIDFilter.SANN:')
msnid <- apply_filter(msnid, filtObj.sa)
show(msnid)
cat('\n')
print('Remove Decoy PSMs:')
msnid <- apply_filter(msnid, "isDecoy == FALSE")
show(msnid)
cat('\n')
print('Remove contaminant PSMs:')
msnid <- apply_filter(msnid, "!grepl('Contaminant',description)")
show(msnid)
cat('\n')

psm.df <- psms(msnid)
psm.df <- add_rownames(psm.df, "Row")

peptide.df <- data.frame(peptides(msnid))
peptide.df <- add_rownames(peptide.df, "Row")

accession.df <- data.frame(accessions(msnid))
accession.df <- add_rownames(accession.df, "Row")


table_dir <- paste(outdir, '/tables',sep='')
dir.create(table_dir, showWarnings=TRUE,recursive=FALSE,mode='0777')
write.table(psm.df, paste(table_dir, '/psms.txt',sep=''),sep='\t', row.names=FALSE)
write.table(peptide.df, paste(table_dir, '/peptides.txt',sep=''),sep='\t',row.names=FALSE)
write.table(accession.df, paste(table_dir, '/accessions.txt',sep=''),sep='\t', row.names=FALSE)

msnset <- as(msnid, "MSnSet")
library("MSnbase")

msnset <- combineFeatures(msnset, fData(msnset)$accession, redundancy.handler="unique", fun="sum",cv=FALSE)
counts.df <- data.frame(exprs(msnset))
counts.df <- add_rownames(counts.df, "Row")
write.table(counts.df, paste(table_dir, '/spectral_counts.txt',sep=''),sep='\t', row.names=FALSE)

sink()
print(head(psms(msnid)))
