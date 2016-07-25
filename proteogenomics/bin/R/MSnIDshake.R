#!/opt/exp_soft/R-3.2.0/bin/Rscript 

library("dplyr")
library("optparse")

option_list = list(
 make_option(c("-i", "--indir"), type="character", default=NULL, 
              help="Directory of mzIdentML files for global FDR control", metavar="character"),
	make_option(c("-v", "--fdr_value"), type="integer", default=1, 
              help="Percentage value to control the FDR", metavar="character"),
	make_option(c("-l","--fdr_level"), type="character", default='accession',
              help="The level to control the FDR ('PSM', 'peptide' or 'accession')", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

fdr_value = opt$fdr_value/100.0   # convert percentage to ratio
fdr_level = opt$fdr_level         

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


#sink(paste(outdir, '/log.txt',sep=''))

print(dirfiles)
library("MSnID")

msnid <- MSnID(".")
msnid <- read_mzIDs(msnid, dirfiles)
msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")

#all_descriptions <- msnid$description

#############
#Unfiltered #
#############

dir.create(paste(outdir,'/qc',sep=''),showWarnings=TRUE,recursive=FALSE,mode='0777')
jpeg(paste(outdir,'/qc/missed_cleavages.jpeg',sep=''))
pepCleav <- unique(psms(msnid)[,c("numMissCleavages", "isDecoy", "peptide")])
pepCleav <- as.data.frame(table(pepCleav[,c("numMissCleavages", "isDecoy")]))
library("ggplot2")
ggplot(pepCleav, aes(x=numMissCleavages, y=Freq, fill=isDecoy)) + geom_bar(stat='identity', position='dodge')  + ggtitle("Number of Missed Cleavages")
dev.off()

jpeg(paste(outdir,'/qc/peptide_lengths.jpeg',sep=''))
msnid$PepLength <- nchar(msnid$peptide) - 4
pepLen <- unique(psms(msnid)[,c("PepLength", "isDecoy", "peptide")])
ggplot(pepLen, aes(x=PepLength, fill=isDecoy)) + geom_histogram(position='dodge', binwidth=3) + ggtitle("Distribution of Peptide Lengths")
dev.off()

#parent ion mass measurement error (ppm)
jpeg(paste(outdir,'/qc/parent_ion_mass_measurement_error.jpeg',sep=''))
ppm <- mass_measurement_error(msnid)
ggplot(as.data.frame(ppm), aes(x=ppm))  + geom_histogram(binwidth=100)
dev.off()

# MSnID package provide a simple correct_peak_selection function that simply adds or subtracts the difference between 13C and 12C to make the error less then 1 Dalton.
jpeg(paste(outdir,'/qc/parent_ion_mass_measurement_error_corrected.jpeg',sep=''))
msnid.fixed <- correct_peak_selection(msnid)
ppm <- mass_measurement_error(msnid.fixed)
ggplot(as.data.frame(ppm), aes(x=ppm)) + geom_histogram(binwidth=0.25)
dev.off()

# recalibrate
jpeg(paste(outdir,'/qc/parent_ion_mass_measurement_error_corrected_recalibrated.jpeg',sep=''))
msnid <- recalibrate(msnid.fixed)
ppm <- mass_measurement_error(msnid)
ggplot(as.data.frame(ppm), aes(x=ppm)) +geom_histogram(binwidth=0.25)
dev.off()

# Shows how decoys distribute with error
jpeg(paste(outdir,'/qc/parent_ion_mass_measurement_error_corrected_recalibrated_stacked.jpeg',sep=''))
temp <- psms(msnid)
temp$error <- ppm
temp$ppm <- temp$error
temp$isDecoy <- msnid$isDecoy
ggplot(temp, aes(x=ppm, fill=isDecoy)) + geom_histogram(position='stack', binwidth=0.25)
dev.off()


# Lets count the number of replicates per peptide
jpeg(paste(outdir,'/qc/peptides_replicate_counts.jpeg',sep=''))
msnid.psms <- psms(msnid)
repcount <- count(msnid.psms, pepSeq, spectrumFile)
repcount <- count(repcount, pepSeq)  # ie. how many replicates per peptide
colnames(repcount)[ncol(repcount)] <- "unf.pepSeq.repcount"
msnid.psms <- merge(msnid.psms, repcount) 

msnid$unf.pepSeq.repcount <- msnid.psms$unf.pepSeq.repcount[match(msnid$pepSeq, msnid.psms$pepSeq)]


repCount <- unique(psms(msnid)[,c("pepSeq", "isDecoy", "unf.pepSeq.repcount")])
repCount <- as.data.frame(table(repCount[,c("unf.pepSeq.repcount", "isDecoy")]))
ggplot(repCount, aes(x=unf.pepSeq.repcount, y=Freq, fill=isDecoy)) + geom_bar(stat='identity', position='dodge') + ggtitle("Distribution of peptide replicate counts")
dev.off()

# PSM condifence (converted to PEP scores) 
msnid$`PeptideShaker PSM confidence` <- as.numeric(msnid$`PeptideShaker PSM confidence`)
msnid$PEP <- 1.0 - msnid$`PeptideShaker PSM confidence`/100.0 

jpeg(paste(outdir,'/qc/psm_PEP.jpeg',sep=''))
params <- psms(msnid)[,c("PEP","isDecoy")]
ggplot(params) + geom_density(aes(x = PEP, color = isDecoy, ..count..))
dev.off()

# PSM mass error - parent
jpeg(paste(outdir,'/qc/corrected_recalibrated_absParentMassErrorPPM.jpeg',sep=''))
msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))
params <- psms(msnid)[,c("absParentMassErrorPPM","isDecoy")]
ggplot(params) + geom_density(aes(x = absParentMassErrorPPM, color = isDecoy, ..count..))
dev.off()

# Parent
jpeg(paste(outdir,'/qc/corr_recal_parent_ion_mass_measurement_targetdecoy.jpeg',sep=''))
dM <- with(psms(msnid),+ (experimentalMassToCharge-calculatedMassToCharge)*chargeState) 
x <- data.frame(dM, isDecoy=msnid$isDecoy)
ggplot(x, aes(x=dM, fill=isDecoy)) +geom_histogram(position='stack', binwidth=0.1)
dev.off()



#msnid <- apply_filter(msnid, "unf.pepSeq.repcount >= 1")
#temp <- transform(psms(msnid), PEPladj = PEP / PepLength)
#msnid$PEPladj <- temp$PEPladj

show(msnid)
msnid$minPepLength <- msnid$PepLength

print('Creating MSnIDFilter object:')

filtObj <- MSnIDFilter(msnid)

filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
filtObj$PEP <- list(comparison="<", threshold=0.01)
filtObj$PepLength <- list(comparison="<=", threshold=35)
filtObj$minPepLength <- list(comparison=">=", threshold=6)
filtObj$numMissCleavages <- list(comparison="<=", threshold=2)
filtObj$numIrregCleavages <- list(comparison="<=", threshold=2)

#filtObj$unf.pepSeq.repcount <- list(comparison=">=", threshold=1)
show(filtObj)
evaluate_filter(msnid, filtObj, level=fdr_level)
cat('\n')

#Protein level FDR filtering
print(paste('Applying ',fdr_level,' level FDR control: ',fdr_value,sep=''))
print('Optimizing MSnIDFilter object with "Grid" method:')
filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max= fdr_value, method="Grid",level=fdr_level, n.iter=1000) 
cat('\n')

print('Optimizing MSnIDFilter.Grid object with "Nelder-Mead" method:')
filtObj.nm <- optimize_filter(filtObj.grid, msnid, fdr.max=fdr_value, method="Nelder-Mead",level=fdr_level, n.iter=1000) 
cat('\n')

print('Optimizing MSnIDFilter.NM object with "Simulated-Annealing" method:')
filtObj.sa <- optimize_filter(filtObj.nm, msnid, fdr.max=fdr_value, method="SANN",level=fdr_level, n.iter=1000) 
cat('\n')

print('Filter msnid object using the MSnIDFilter.SANN:')
msnid <- apply_filter(msnid, filtObj.sa)


show(filtObj.sa)
show(msnid)

print('...DONE APPLYING FDR STRINGENCY CRITERIA...')

# Done applying FDR filtering criteria

jpeg(paste(outdir,'/filtered_PSM_PEP.jpeg',sep=''))
params <- psms(msnid)[,c("PEP","isDecoy")]
ggplot(params) + geom_density(aes(x = PEP, color = isDecoy, ..count..))
dev.off()


decoy <- apply_filter(msnid, "isDecoy == TRUE")
contaminants <- apply_filter(msnid, "grepl('CONTAMINANT',description)")

msnid <- apply_filter(msnid, "isDecoy == FALSE")
#msnid <- apply_filter(msnid, "!grepl('CONTAMINANT',description)")

contaminant.df <- psms(contaminants)
contaminant.df <- add_rownames(contaminant.df, "Row")

decoy.df <- psms(decoy)
decoy.df <- add_rownames(decoy.df, "Row")

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

write.table(contaminant.df, paste(table_dir, '/contaminants.txt',sep=''),sep='\t', row.names=FALSE)
write.table(decoy.df, paste(table_dir, '/decoys.txt',sep=''),sep='\t', row.names=FALSE)


msnset <- as(msnid, "MSnSet")
library("MSnbase")

msnset <- combineFeatures(msnset, fData(msnset)$accession, redundancy.handler="unique", fun="sum",cv=FALSE)
counts.df <- data.frame(exprs(msnset))
counts.df <- add_rownames(counts.df, "Row")
write.table(counts.df, paste(table_dir, '/spectral_counts.txt',sep=''),sep='\t', row.names=FALSE)





