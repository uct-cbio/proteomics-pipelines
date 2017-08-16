#!/usr/bin/env bash

#PBS -N HumanFastaExport
#PBS -q UCTlong
#PBS -l nodes=1:ppn=10:series600
#PBS -M matthys.potgieter@gmail.com
#PBS -m ae

# NB, PLEASE READ THIS!
# There is a 2:1 correspondence between RAM and cores on the 600 series.
# You need to know how much RAM your job will consume before submitting it.
# Please set the ppn value above to be 1/2 the GB of RAM required.  For
# example a job needing 10GB of RAM should have ppn=5

# Please leave the hostname command here for troubleshooting purposes.
hostname

set -e

cd /scratch/DB/bio/metaproteomics/taxonomic_divisions_uniprot

#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/reldate.txt
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.dat.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_bacteria.dat.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_fungi.dat.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_human.dat.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_viruses.dat.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_archaea.dat.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_bacteria.dat.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_fungi.dat.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_human.dat.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_viruses.dat.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/README

#gunzip *.gz
varsplic.pl -input uniprot_sprot_human.dat     -check_vsps -crosscheck -error varsplic_sprot_human.err     -fasta uniprot_sprot_human_varsplic.fasta   -which full
varsplic.pl -input uniprot_trembl_human.dat    -check_vsps -crosscheck -error varsplic_trembl_human.err    -fasta uniprot_trembl_human_varsplic.fasta   -which full

#varsplic.pl -input uniprot_sprot_archaea.dat   -check_vsps -crosscheck -error varsplic_sprot_archaea.err   -fasta uniprot_sprot_archaea_varsplic.fasta   -which full
#varsplic.pl -input uniprot_sprot_bacteria.dat  -check_vsps -crosscheck -error varsplic_sprot_bacteria.err  -fasta uniprot_sprot_bacteria_varsplic.fasta  -which full
#varsplic.pl -input uniprot_sprot_fungi.dat     -check_vsps -crosscheck -error varsplic_sprot_fungi.err     -fasta uniprot_sprot_fungi_varsplic.fasta     -which full
#varsplic.pl -input uniprot_sprot_viruses.dat   -check_vsps -crosscheck -error varsplic_sprot_viruses.err   -fasta uniprot_sprot_viruses_varsplic.fasta   -which full

#varsplic.pl -input uniprot_trembl_archaea.dat  -check_vsps -crosscheck -error varsplic_trembl_archaea.err  -fasta uniprot_trembl_archaea_varsplic.fasta   -which full
#varsplic.pl -input uniprot_trembl_bacteria.dat -check_vsps -crosscheck -error varsplic_trembl_bacteria.err -fasta uniprot_trembl_bacteria_varsplic.fasta  -which full
#varsplic.pl -input uniprot_trembl_fungi.dat    -check_vsps -crosscheck -error varsplic_trembl_fungi.err    -fasta uniprot_trembl_fungi_varsplic.fasta     -which full
#varsplic.pl -input uniprot_trembl_viruses.dat  -check_vsps -crosscheck -error varsplic_trembl_viruses.err  -fasta uniprot_trembl_viruses_varsplic.fasta   -which full
