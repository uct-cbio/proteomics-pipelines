#!/usr/bin/env bash
set -e
JVM_ARGS="-d64 -Xms1024M -Xmx15360M -server"
configfile=$1
eval "$(egrep '^#|^[^ ]*=[^;&]*'  "$configfile")"
########
rm -rf ${output_folder}
mkdir ${output_folder}

cp $configfile $output_folder/config.cfg

mzid_folder=${output_folder}"/mzIdentMLs"
mkdir ${mzid_folder}

fasta_folder=${output_folder}"/fasta"
mkdir ${fasta_folder}

cp ${target_fasta} ${fasta_folder}
cp ${contaminant_fasta} ${fasta_folder}

contaminant=$(basename $contaminant_fasta)
target=$(basename $target_fasta)

contaminant=${contaminant%.fasta}
target=${target%.fasta}

fasta=${fasta_folder}/${target}_${contaminant}.fasta
cat ${target_fasta} <(annotateContaminantUniprotFasta.py ${contaminant_fasta}) > ${fasta}

unvalidated_folder=${output_folder}"/unvalidated"

mkdir ${unvalidated_folder}

temp_folder=${output_folder}"/temp"
log_folder=${output_folder}"/log"

mkdir ${temp_folder}
mkdir ${log_folder}

cp -R ${ps_folder} ${output_folder}/ps
cp -R ${sg_folder} ${output_folder}/sg
wait

# Set paths for logs and temporary files
java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.PathSettingsCLI -temp_folder ${temp_folder} -log ${log_folder} 
wait

# create target-decoy fasta in the output directory
java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar \
eu.isas.searchgui.cmd.FastaCLI \
-in ${fasta} -decoy 
wait

# Create search parameters
java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -out ${output_folder}/identification.par -db ${fasta%.fasta}_concatenated_target_decoy.fasta -prec_tol ${prec_tol} -prec_ppm ${prec_ppm} -frag_tol ${frag_tol} -frag_ppm ${frag_ppm} -enzyme "${enzyme}" -fixed_mods "${fixed_mods}" -variable_mods "${variable_mods}" -min_charge ${min_charge} -max_charge ${max_charge} -mc ${mc} -fi ${fi} -ri ${ri} -psm_fdr ${psm_fdr} -peptide_fdr ${peptide_fdr} -protein_fdr ${protein_fdr} -myrimatch_min_pep_length ${myrimatch_min_pep_length} -myrimatch_max_pep_length ${myrimatch_max_pep_length} -msgf_instrument ${msgf_instrument} -msgf_min_pep_length ${msgf_min_pep_length} -msgf_max_pep_length ${msgf_max_pep_length} -tide_min_pep_length ${tide_min_pep_length} -tide_max_pep_length ${tide_max_pep_length} -import_peptide_length_min ${import_peptide_length_min} -import_peptide_length_max ${import_peptide_length_max} -annotation_level ${annotation_level} 
wait

