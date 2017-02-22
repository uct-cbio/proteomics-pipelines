#!/usr/bin/env bash

config=$1
outpath=$2
kegg_id=$3 #'mtu'
python2ve=$4 # $HOME/ve273   # interproscan needs python 2 (Create one with virtualenv)

# The rest is generic

# Create the output directory
#rm -rf $outpath && mkdir $outpath

# Check the config script for any errors
#mq_proteogenomics_validate_config.py $config $outpath

# Proteogenomics analysis
if [ ! -d $outpath/uniprot ] ; then
    mq_uniprot.py $config $outpath || rm -rf $outpath/uniprot
fi

if [ ! -d $outpath/strains ] ; then
    mq_genome_to_peptide.py $config $outpath #|| rm -rf $outpath/strains
fi

#uniprot_peptide2db.py $config $outpath

if [ ! -d $outpath/mapping ] ; then
    mq_peptide_to_referencedb.py $config $outpath || rm -rf $outpath/mapping
fi

if [ ! -d $outpath/fasta ] ; then
    mq_fasta_export.py $config $outpath || rm -rf $outpath/fasta
fi

if [ ! -f $outpath/fasta/id_mapping.p ] ; then
    ips_fasta.py $outpath/fasta/combined_translated.fasta $outpath/fasta || rm -rf $outpath/fasta/id_mapping.p
fi

if [ ! -d $outpath/blast ] ; then
    mkdir $outpath/blast
    mq_blast2ref.py $config $outpath  || rm -rf $outpath/blast 
fi

#if [ ! -f $outpath/fasta/nr_translated_pg_orfs.fasta.gff3 ] ; then
#    ips.sh $outpath/fasta/nr_translated_pg_orfs.fasta $outpath/fasta $python2ve
#fi

if [ ! -f $outpath/combined.csv ] ; then
    mq_peptide_to_protein.py $config $outpath || rm -rf $outpath/combined.csv
fi

if [ ! -d $outpath/tables ] ;  then      
    mq_export_tables.py $config $outpath || rm -rf $outpath/tables
fi

#############################
# Identification Statistics #
#############################

if [ ! -d $outpath/stats ] ; then
    mq_basestats.py $config $outpath || rm -rf $outpath/stats
fi

#################
# GBROWSE DATA  #
#################
mq_peptide_features.py $config $outpath
#mq_domain_features.py $config $outpath
#mq_contig_heatmaps.py $config $outpath
#mq_wiggle_features.py $config $outpath

##################
# IPS Enrichment #
##################
#ips_gsea.py $outpath 
#mq_annotate.py $outpath
#mq_genesets.R --outdir $outpath --keggid $kegg_id

###################################
# Differential abundance analysis #
###################################
#mq_experimental_design.py $config $outpath
#rm -rf ${outpath}/diff
#mq_differential_abundance.R -d ${outpath}/experimental_design.R -p ${outpath}/combined.csv -o ${outpath}/diff






