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


# Proteogenomics analysis
# check the ncrna_class tag in GFF3
if [ ! -d $outdir/ena ] ; then
    mq_ena.py $config $outdir || ( rm -rf $outdir/ena ; exit 1 )
fi



if [ ! -d $outdir/strains ] ; then
    mq_genome_to_peptide.py $config $outdir || ( rm -rf $outdir/strains ; exit 1 )
fi


#stub_mq_consensus_blast.py $config $outdir

#uniprot_peptide2db.py $config $outdir

outfile=$outdir/mapping 
if [ ! -d $outfile ] ; then
    mq_peptide_to_referencedb.py $config $outdir || ( rm -rf $outfile ; exit 1 )
fi


#########
# BLAST #
#########

if [ ! -d $outdir/blast ] ; then
    mkdir $outdir/blast
fi

################################################ 
# BLAST or protein group FASTA to SIX FRAME DB #
################################################
# speed up by using an nr database for query fasta
# speed up by using an nr database for target fasta

# map the search sequence to the stop-to-stop orf sequence, uses the top match per orf
if [ ! -d $outdir/blast/groups2orfs ] ; then
    mkdir $outdir/blast/groups2orfs
    mq_blast_groups2orfs.py $config $outdir  || ( rm -rf $outdir/blast/groups2orfs ; exit 1 )
fi

# map the orfs to the proteomes, top match per orf exported
if [ ! -d $outdir/blast/orfs2proteins ] ; then
    mkdir $outdir/blast/orfs2proteins
    mq_blast_orfs2refproteome.py $config $outdir  || ( rm -rf $outdir/blast/orfs2proteins ; exit 1 )
fi


if [ ! -d $outdir/blast/orfs2genome ] ; then
    mkdir $outdir/blast/orfs2genome
    mq_blast_orfs2refgenome.py $config $outdir  || ( rm -rf $outdir/blast/orfs2genome ; exit 1 )
fi

# this output is not used for now, peptides are all mapped to the same reference, to facilitate gbrowse on the same reference sequence if it is needed
#if [ ! -d $outdir/blast/peptides2genome ] ; then
#    mkdir $outdir/blast/peptides2genome
#    mq_blast_peptides2refgenome.py $config $outdir  #|| rm -rf $outdir/blast/peptides2genome
#fi

if [ ! -d $outdir/blast/peptides2orfs ] ; then
    mkdir $outdir/blast/peptides2orfs
    mq_blast_peptides2reforfs.py $config $outdir  || ( rm -rf $outdir/blast/peptides2orfs ; exit 1 )
fi

##################
# Operon Mapping #
##################

#if [ ! -f $outdir/mapping/operons.json ] ; then
#    door2_operon_map.py $config $outdir || ( rm -rf $outdir/mapping/operons.json ; exit 1 )
#fi

##################
# Combined table #
##################
if [ ! -f $outdir/config.yml ] ; then
    echo $outdir/config.yml
    mqmetaproteomics.py $config
fi
#if [ ! -f $outdir/combined.csv ] ; then
#ls $outdir/*combined.csv > /dev/null 2>&1  || mq_peptide_to_protein.py $config $outdir || ( rm -rf $outdir/*combined.csv ; exit 1 )

#############################
# Identification Statistics #
#############################
if [ ! -d $outdir/stats ] ; then
    mq_basestats.py $config $outdir || rm -rf $outdir/stats
    mq_annotations.py $config $outdir || rm -rf $outdir/stats
fi
#mq_annotations_ko.py $config $outdir

mq_annotations_pathview.R --config $config --outdir $outdir

echo "Annotations exported"
exit 0

#if [ ! -d $outdir/tables2 ] ;  then      
#    mq_export_tables.py $config $outdir || rm -rf $outdir/tables
#fi

#exit 0

if [ ! -d $outdir/jbrowse ] ; then
    mkdir $outdir/jbrowse
    mq_strains2ref.py $config $outdir      || ( rm -rf $outdir/jbrowse ; exit 1 )
    # not needed  #mq_strains2ref_peptides.py $config $outdir  || ( rm -rf $outdir/jbrowse ; exit 1 )
fi

mq_peptide_features.py $config $outdir
#jbrowse_iprscan.py $config $outdir

mq_jbrowse_upload_script.py $config $outdir || ( rm -rf $outdir/jbrowse_ ; exit 1 )

cd $outdir/jbrowse/local && echo $hostname && jbrowse admin-server

#mq_contig_heatmaps.py $config $outdir
#mq_wiggle_features.py $config $outdir








###################### 
# PEP score analysis #
######################

#mq_PEP_MSMS.py $config $outdir
#mq_differential_abundance.R -d ${outdir}/experimental_design.R -p ${outdir}/combined.csv -o ${outdir}/diff # || rm -rf ${outdir}/diff




