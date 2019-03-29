#!/usr/bin/env bash

set -e

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
    mq_genome_to_peptide.py $config $outpath || rm -rf $outpath/strains
fi

#mq_consensus_blast.py $config $outpath

#uniprot_peptide2db.py $config $outpath

################
# Interproscan #
################

if [ ! -d $outpath/mapping ] ; then
    mq_peptide_to_referencedb.py $config $outpath || rm -rf $outpath/mapping
fi


if [ ! -d $outpath/fasta ] ; then
    mq_fasta_export.py $config $outpath || rm -rf $outpath/fasta
fi


if [ ! -f $outpath/fasta/id_mapping.json ] ; then
    ips_fasta.py $outpath/fasta/combined_translated.fasta $outpath/fasta || rm -rf $outpath/fasta/id_mapping.json
fi


if [ ! -f $outpath/fasta/nr_translated_pg_orfs.fasta.gff3 ] ; then
    ips.sh $outpath/fasta/nr_translated_pg_orfs.fasta $outpath/fasta $python2ve || rm -rf $outpath/fasta/nr_translated_pg_orfs.fasta.gff3
fi


#########
# BLAST #
#########

if [ ! -d $outpath/blast ] ; then
    mkdir $outpath/blast
fi


if [ ! -d $outpath/blast/orfs2proteins ] ; then
    mkdir $outpath/blast/orfs2proteins
    mq_blast_orfs2refproteome.py $config $outpath  || rm -rf $outpath/blast/orfs2proteins
fi



if [ ! -d $outpath/blast/orfs2genome ] ; then
    mkdir $outpath/blast/orfs2genome
    mq_blast_orfs2refgenome.py $config $outpath  || rm -rf $outpath/blast/orfs2genome
fi

#if [ ! -d $outpath/blast/peptides2genome ] ; then
#    mkdir $outpath/blast/peptides2genome
#    mq_blast_peptides2refgenome.py $config $outpath  || rm -rf $outpath/blast/peptides2genome
#fi

if [ ! -d $outpath/blast/peptides2orfs ] ; then
    mkdir $outpath/blast/peptides2orfs
    mq_blast_peptides2reforfs.py $config $outpath  || rm -rf $outpath/blast/peptides2orfs
fi

##################
# Operon Mapping #
##################

if [ ! -f $outpath/mapping/operons.json ] ; then
    door2_operon_map.py $config $outpath || rm -rf $outpath/mapping/operons.json
fi

##################
# Combined table #
##################

if [ ! -f $outpath/combined.csv ] ; then
    mq_peptide_to_protein.py $config $outpath || rm -rf $outpath/combined.csv ; exit 1
fi

###################################
# Differential abundance analysis #
###################################

if [ ! -f $outpath/experimental_design.R ] ; then
    mq_experimental_design.py $config $outpath || rm -rf $outpath/experimental_design.R
fi

if [ ! -d $outpath/msnbase ] ; then
    mq_normalize_intensity.R -q 'iBAQ.'  -p ${outpath}/combined.csv -o ${outpath} || rm -rf ${outpath}/msnbase
fi


if [ ! -d $outpath/diff ] ; then
    mq_diff.R -d ${outpath}/experimental_design.R -f ${outpath}/msnbase/normalized.csv -o ${outpath}/diff || rm -rf ${outpath}/diff
    #mq_differential_abundance.R -d ${outpath}/experimental_design.R -p ${outpath}/combined.csv -o ${outpath}/diff || rm -rf ${outpath}/diff
fi
##################
# IPS Enrichment #
##################
if [ ! -d $outpath/gsea ] ; then
    ips_gsea.py $outpath  $kegg_id #&& mq_annotate.py $outpath && mq_genesets.R --outdir $outpath/gsea --keggid $kegg_id #|| rm -rf $outpath/gsea ; exit 1
    gage.R --outdir $outpath/gsea --indir $outpath/gsea --keggid $kegg_id --design ${outpath}/experimental_design.R --table ${outpath}/msnbase/normalized.csv --genecol GeneID --kocol GeneID --pval 0.05 || rm -rf ${outpath}/gsea
fi

exit

###########
# JBROWSE #
###########

if [ ! -d $outpath/jbrowse ] ; then
    mkdir $outpath/jbrowse
    mq_strains2ref_peptides.py $config $outpath  || rm -rf $outpath/jbrowse
    mq_strains2ref_orfs.py $config $outpath      || rm -rf $outpath/jbrowse
    mq_jbrowse_upload_script.py $config $outpath || rm -rf $outpath/jbrowse
fi


#mq_strains2ref_orfs.py $config $outpath      

#mq_peptide_features.py $config $outpath
#mq_domain_features.py $config $outpath
#mq_contig_heatmaps.py $config $outpath
#mq_wiggle_features.py $config $outpath





#if [ ! -d $outpath/tables ] ;  then      
#    mq_export_tables.py $config $outpath || rm -rf $outpath/tables
#fi

#############################
# Identification Statistics #
#############################
#if [ ! -d $outpath/stats ] ; then
#    mq_basestats.py $config $outpath || rm -rf $outpath/stats
#fi


###################### 
# PEP score analysis #
######################

#mq_PEP_MSMS.py $config $outpath
#mq_differential_abundance.R -d ${outpath}/experimental_design.R -p ${outpath}/combined.csv -o ${outpath}/diff # || rm -rf ${outpath}/diff




