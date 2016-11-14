#!/usr/bin/env bash

set -e

cd /scratch/DB/bio/metaproteomics/taxonomic_divisions_uniprot

wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/reldate.txt
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_bacteria.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_fungi.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_human.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_viruses.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_archaea.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_bacteria.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_fungi.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_human.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_viruses.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/README

gunzip *

perl varsplic.pl -input uniprot_sprot_archaea.dat.gz  -check_vsps -crosscheck -error varsplic.err -fasta uniprot_sprot_archaea_varsplic.fasta -which full
perl varsplic.pl -input uniprot_sprot_bacteria.dat.gz -check_vsps -crosscheck -error varsplic.err -fasta uniprot_sprot_bacteria_varsplic.fasta -which full
perl varsplic.pl -input uniprot_sprot_fungi.dat.gz    -check_vsps -crosscheck -error varsplic.err -fasta uniprot_sprot_fungi_varsplic.fasta -which full
perl varsplic.pl -input uniprot_sprot_viruses.dat.gz  -check_vsps -crosscheck -error varsplic.err -fasta uniprot_sprot_viruses_varsplic.fasta -which full

perl varsplic.pl -input uniprot_trembl_archaea.dat.gz  -check_vsps -crosscheck -error varsplic.err -fasta uniprot_trembl_archaea_varsplic.fasta -which full
perl varsplic.pl -input uniprot_trembl_bacteria.dat.gz -check_vsps -crosscheck -error varsplic.err -fasta uniprot_trembl_bacteria_varsplic.fasta -which full
perl varsplic.pl -input uniprot_trembl_fungi.dat.gz    -check_vsps -crosscheck -error varsplic.err -fasta uniprot_trembl_fungi_varsplic.fasta -which full
perl varsplic.pl -input uniprot_trembl_viruses.dat.gz  -check_vsps -crosscheck -error varsplic.err -fasta uniprot_trembl_viruses_varsplic.fasta -which full
