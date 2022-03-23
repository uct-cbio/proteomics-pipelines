#!/usr/bin/env bash

set -e
set -a 

config=$1
: ${1?"Please provide the path to the yaml config file"}



outdir=$(extract_variable.py $config outdir) || ( echo "no outdir defined"  && exit 1)

#kegg_id=$3 #'mtu'
#: ${3?"Please provide the kegg id for this species eg : mtu"}

normalize=$(extract_variable.py $config normalize) || ( echo "no 'normalize' method defined eg. 'quantiles'"  && exit 1)

impute=$(extract_variable.py $config impute) || ( echo "no 'impute' method defined eg. 'knn'"  && exit 1)

# Proteogenomics analysis
if [ ! -d $outdir/uniprot ] ; then
    mq_uniprot.py $config $outdir || ( rm -rf $outdir/uniprot ; exit 1 )
fi

mqmetaproteomics.py $config

if [ ! -d $outdir/strains ] ; then
    mq_genome_to_peptide.py $config $outdir || ( rm -rf $outdir/strains ; exit 1 )
fi


#mq_consensus_blast.py $config $outdir


#uniprot_peptide2db.py $config $outdir

exit 0
################
# Interproscan #
################

outfile=$outdir/mapping 
if [ ! -d $outfile ] ; then
    mq_peptide_to_referencedb.py $config $outdir || ( rm -rf $outfile ; exit 1 )
fi


if [ ! -d $outdir/fasta ] ; then
    mq_fasta_export.py $config $outdir || ( rm -rf $outdir/fasta ; exit 1 )
fi


if [ ! -f $outdir/fasta/id_mapping.json ] ; then
    ips_fasta.py $outdir/fasta/combined_translated.fasta $outdir/fasta || ( rm -rf $outdir/fasta/id_mapping.json ; exit 1 )
fi

# Run InterProScan
outfile=$outdir/fasta/nr_translated_pg_orfs.fasta.gff3 
if [ ! -f $outfile ] ; then
    ips.sh $outdir/fasta/nr_translated_pg_orfs.fasta $outdir/fasta || ( rm -rf $outfile ; echo failed  ; exit 1 )
fi


#########
# BLAST #
#########

if [ ! -d $outdir/blast ] ; then
    mkdir $outdir/blast
fi


if [ ! -d $outdir/blast/orfs2proteins ] ; then
    mkdir $outdir/blast/orfs2proteins
    mq_blast_orfs2refproteome.py $config $outdir  || ( rm -rf $outdir/blast/orfs2proteins ; exit 1 )
fi



if [ ! -d $outdir/blast/orfs2genome ] ; then
    mkdir $outdir/blast/orfs2genome
    mq_blast_orfs2refgenome.py $config $outdir  || ( rm -rf $outdir/blast/orfs2genome ; exit 1 )
fi

#if [ ! -d $outdir/blast/peptides2genome ] ; then
#    mkdir $outdir/blast/peptides2genome
#    mq_blast_peptides2refgenome.py $config $outdir  || rm -rf $outdir/blast/peptides2genome
#fi

if [ ! -d $outdir/blast/peptides2orfs ] ; then
    mkdir $outdir/blast/peptides2orfs
    mq_blast_peptides2reforfs.py $config $outdir  || ( rm -rf $outdir/blast/peptides2orfs ; exit 1 )
fi

##################
# Operon Mapping #
##################

if [ ! -f $outdir/mapping/operons.json ] ; then
    door2_operon_map.py $config $outdir || ( rm -rf $outdir/mapping/operons.json ; exit 1 )
fi

##################
# Combined table #
##################

if [ ! -f $outdir/combined.csv ] ; then
    mq_peptide_to_protein.py $config $outdir || ( rm -rf $outdir/combined.csv ; exit 1 )
fi

###################################
# Differential abundance analysis #
###################################

if [ ! -f $outdir/experimental_design.R ] ; then
    mq_experimental_design.py $config $outdir || ( rm -rf $outdir/experimental_design.R ; exit 1 )
fi

if [ ! -d $outdir/msnbase ] ; then
    mq_normalize_intensity.R -q 'iBAQ.'  -p ${outdir}/combined.csv -o ${outdir} -n ${normalize} -i ${impute} || ( rm -rf ${outdir}/msnbase ; exit 1 )
fi

if [ ! -d $outdir/diff ] ; then
    mq_diff.R -d ${outdir}/experimental_design.R -f ${outdir}/msnbase/normalized.csv -o ${outdir}/diff || ( rm -rf ${outdir}/diff ; exit 1 )
    #mq_differential_abundance.R -d ${outdir}/experimental_design.R -p ${outdir}/combined.csv -o ${outdir}/diff || rm -rf ${outdir}/diff
fi

##################
# IPS Enrichment #
##################
if [ ! -d $outdir/gsea ] ; then
    ips_gsea.py $outdir  $kegg_id && gage.R --outdir $outdir/gsea --indir $outdir/gsea --keggid $kegg_id --design ${outdir}/experimental_design.R --table ${outdir}/msnbase/normalized.csv --genecol GeneID --kocol GeneID --pval 0.05 || ( rm -rf ${outdir}/gsea ; exit 1 )
fi
echo 'Done'
exit
###########
# JBROWSE #
###########

if [ ! -d $outdir/jbrowse ] ; then
    mkdir $outdir/jbrowse
    mq_strains2ref_peptides.py $config $outdir  || ( rm -rf $outdir/jbrowse ; exit 1 )
    mq_strains2ref_orfs.py $config $outdir      || ( rm -rf $outdir/jbrowse ; exit 1 )
    mq_jbrowse_upload_script.py $config $outdir || ( rm -rf $outdir/jbrowse ; exit 1 )
fi


#mq_strains2ref_orfs.py $config $outdir      

#mq_peptide_features.py $config $outdir
#mq_domain_features.py $config $outdir
#mq_contig_heatmaps.py $config $outdir
#mq_wiggle_features.py $config $outdir





#if [ ! -d $outdir/tables ] ;  then      
#    mq_export_tables.py $config $outdir || rm -rf $outdir/tables
#fi

#############################
# Identification Statistics #
#############################
#if [ ! -d $outdir/stats ] ; then
#    mq_basestats.py $config $outdir || rm -rf $outdir/stats
#fi


###################### 
# PEP score analysis #
######################

#mq_PEP_MSMS.py $config $outdir
#mq_differential_abundance.R -d ${outdir}/experimental_design.R -p ${outdir}/combined.csv -o ${outdir}/diff # || rm -rf ${outdir}/diff




