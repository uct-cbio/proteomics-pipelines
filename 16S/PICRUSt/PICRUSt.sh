#!/bin/bash
#PBS -N picrust
#PBS -S /bin/bash
#PBS -q UCTlong
#PBS -l nodes=1:series600:ppn=4
#PBS -V
#PBS -d /specify/directory

source /home/kviljoen/activate_qiime.sh #change as appropriate

#NB if having issues with calling the wrong python disable the above line
inDir=/specify/input_dir
outDir=$inDir #specify output directory

/opt/exp_soft/python-2.7.3/bin/normalize_by_copy_number.py -i $inDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline.biom -o $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline.CN_norm.biom

/opt/exp_soft/python-2.7.3/bin/predict_metagenomes.py -i $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline.CN_norm.biom -o $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_predicted_metagenome.biom --accuracy_metrics $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_predicted_metagenome.NSTIs.txt 

#with_confidence returns metagenomic 95% confidence intervals for each gene prediction
#the --accuracy_metrics flags returns NSTI scores (The NSTI scores that you get back are the branch length on the Greengenes tree separating taxa to be predicted from the nearest sequenced genome. This set of these distances is averaged, weighting by the abundance of each organism)
#Note: used GG13.8 to map our IDs to greengenes but PICRUSt uses 13.5. This appears to be fine however as stated by one of the PICRUSt authors Daniel Mcdonald here https://groups.google.com/forum/#!topic/picrust-users/LTxmapB1xiA

/opt/exp_soft/python-2.7.3/bin/categorize_by_function.py -i $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_predicted_metagenome.biom -l 3 -c KEGG_Pathways -o $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_categorize_by_function_level_3.biom
/opt/exp_soft/python-2.7.3/bin/categorize_by_function.py -i $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_predicted_metagenome.biom -l 2 -c KEGG_Pathways -o $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_categorize_by_function_level_2.biom
/opt/exp_soft/python-2.7.3/bin/categorize_by_function.py -i $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_predicted_metagenome.biom -l 1 -c KEGG_Pathways -o $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_categorize_by_function_level_1.biom

#convert .biom files to .txt files for use in R
biom convert -i $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_predicted_metagenome.biom -o $outDir/v2_GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_predicted_metagenome.txt -b #b specifies biom --> txt
biom convert -i $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_categorize_by_function_level_3.biom -o $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_categorize_by_function_level_3.txt -b
biom convert -i $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_categorize_by_function_level_2.biom -o $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_categorize_by_function_level_2.txt -b
biom convert -i $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_categorize_by_function_level_1.biom -o $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_categorize_by_function_level_1.txt -b

#STAMP software can also use .txt files the headers just need minor modifications I think (see exisiting .txt files previously used for STAMP) 

#lets see which OTUs underly the predicted KEGG functions using metagenome_contributions.py
/opt/exp_soft/python-2.7.3/bin/metagenome_contributions.py -i $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline.CN_norm.biom -o $outDir/GG_13_8_closed_reference_OTU_table_from_de_novo_pipeline_metagenome_contributions.txt -t ko







