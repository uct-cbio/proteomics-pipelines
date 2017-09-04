#!/bin/bash
#PBS -N uparse_closed
#PBS -S /bin/bash
#PBS -q UCTlong
#PBS -l nodes=1:series600:ppn=4
#PBS -V
#PBS -M email@adress
#PBS -d /specify/directory


source /home/kviljoen/activate_qiime.sh #change as appropriate
inDir=/input_directory/uparse_downstream/ #specify input directory
outDir=/output_directory/ #specify output directory

usearch9 -usearch_global $inDir/otus_repsetOUT.fa -db /scratch/DB/bio/qiime/greengenes/gg_13_8_otus/rep_set/97_otus.fasta  -id 0.97 -strand plus -uc $outDir/de_novo_repset_to_GG_13_8_map.uc 
#NOTE1: now download the resulting .uc file which maps de novo IDs to GG IDs and OTU rownames with GG IDs (for those that map) in R
#NOTE2: next upload the .txt filtered otu table with GG IDs

#convert uploaded otu .txt table to .biom format - the output will be used for PICRUSt
biom convert -i $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline.txt -o $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline.biom --table-type="otu table"
#

#summarize .biom otu table
biom summarize-table -i $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline.biom -o $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline.biom.summary
