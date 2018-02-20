#!/bin/sh

sf_file=$1
read -p "Enter EBI genome accession(s) for alignment: " genome_accession

#genome_accession=$2
cp $sf_file ${sf_file##*/}
sf_file=${sf_file##*/}

genome_accession=`echo "${genome_accession}"|sed -e 's/ /_/g'`
genome ${genome_accession}
genome_name=${genome_accession}.fasta
GN=${genome_name%.*}
SF=${sf_file%.*}
bwa index -p ${GN} -a is $genome_name
bwa bwasw ${GN} ${sf_file}  > ${SF}.sam
samtools faidx $genome_name
samtools import ${genome_name}.fai ${SF}.sam ${SF}.bam
samtools sort ${SF}.bam ${SF}.sort
samtools index ${SF}.sort.bam
samtools flagstat ${SF}.sort.bam|tee flagstat_${GN}.txt





