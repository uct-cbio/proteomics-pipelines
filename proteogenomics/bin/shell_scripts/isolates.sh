#!/bin/bash

#read -p "Enter EBI genome accession for primary alignment: " accession_1
#read -p "Enter EBI genome accession for seconday alignment of non-aligend reads: " accession_2

accession_1='AL123456'    #h37rv
accession_2='AE000516'    #cdc1551
accession_3='CP010873'    #beijing-like   (no reference proteome)

rp_1='UP000001584'
rp_2='UP000001020'

contig_file=$1

cp $contig_file ${contig_file##*/}

contig_file=${contig_file##*/}

sf <(cat $contig_file)>six_frame.fasta

##################  H37rv  ####################################
#aln six_frame.fasta $accession_1
sf_file=six_frame.fasta

genome ${accession_1}
rp ${rp_1} > ${rp_1}.fasta

#genome_accession=`echo "${genome_accession}"|sed -e 's/ /_/g'` uncomment this to align to more than one genome at a time
genome_name=${accession_1}.fasta
GN=${genome_name%.*}
SF=${sf_file%.*}
bwa index -p ${GN} -a is $genome_name
bwa bwasw ${GN} ${sf_file}  > ${SF}.sam
samtools faidx $genome_name
samtools import ${genome_name}.fai ${SF}.sam ${SF}.bam
samtools sort ${SF}.bam ${SF}.sort
samtools index ${SF}.sort.bam
samtools flagstat ${SF}.sort.bam|tee flagstat_${GN}.txt
samtools view -F4 ${SF}.sort.bam|sam2fasta > mapped_${GN}.fasta

cat mapped_${GN}.fasta|nuc2prot ${accession_1}.fasta ${rp_1}.fasta > mapped_${GN}_proteins.fasta
samtools view -f4 ${SF}.sort.bam|sam2fasta > unmapped_${GN}.fasta

rm -rf ${SF}*
rm -rf ${GN}*

########################  CDC 1551  ########################

sf_file=unmapped_${GN}.fasta

genome ${accession_2}
rp ${rp_2} > ${rp_2}.fasta

genome_name=${accession_2}.fasta
GN=${genome_name%.*}
SF=${sf_file%.*}
bwa index -p ${GN} -a is $genome_name
bwa bwasw ${GN} ${sf_file}  > ${SF}.sam
samtools faidx $genome_name
samtools import ${genome_name}.fai ${SF}.sam ${SF}.bam
samtools sort ${SF}.bam ${SF}.sort
samtools index ${SF}.sort.bam
samtools flagstat ${SF}.sort.bam|tee flagstat_${GN}.txt
samtools view -F4 ${SF}.sort.bam|sam2fasta > mapped_${GN}.fasta

cat mapped_${GN}.fasta|nuc2prot ${accession_2}.fasta ${rp_2}.fasta > mapped_${GN}_proteins.fasta
samtools view -f4 ${SF}.sort.bam|sam2fasta > unmapped_${GN}.fasta

rm -rf ${SF}*
rm -rf ${GN}*

###  BEIJING  ##############

sf_file=unmapped_${GN}.fasta
genome ${accession_3}
genome_name=${accession_3}.fasta
GN=${genome_name%.*}
SF=${sf_file%.*}

bwa index -p ${accession_3} -a is $genome_name
bwa bwasw ${accession_3} ${sf_file}  > ${SF}.sam
samtools faidx $genome_name
samtools import ${genome_name}.fai ${SF}.sam ${SF}.bam
samtools sort ${SF}.bam ${SF}.sort
samtools index ${SF}.sort.bam
samtools flagstat ${SF}.sort.bam|tee flagstat_${accession_3}.txt
samtools view -F4 ${SF}.sort.bam|sam2fasta > mapped_${accession_3}.fasta

samtools view -f4 ${SF}.sort.bam|sam2fasta > unmapped_${accession_3}.fasta

cat mapped_${accession_3}.fasta|nuc2prot > mapped_${accession_3}_proteins.fasta
cat unmapped_${accession_3}.fasta|nuc2prot > unmapped_${accession_3}_proteins.fasta

rm -rf ${SF}*
rm -rf ${GN}*
#rp_header='uniprot_protein_sequence'

pg_header='proteogenomics_six_frame_sequence'
#cat ${rp_1}.fasta ${rp_2}.fasta |nr_fasta ${rp_header} > rp_db.fasta

cat mapped_${accession_1}_proteins.fasta mapped_${accession_2}_proteins.fasta mapped_${accession_3}_proteins.fasta  unmapped_${accession_3}_proteins.fasta | nr_fasta ${pg_header} > pg_db.fasta

cat pg_db.fasta |mstart > pg_db_m.fasta

map_fasta pg_db_m.fasta pg_db.fasta ${rp_1}.fasta ${rp_2}.fasta








