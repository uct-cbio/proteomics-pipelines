#!/usr/bin/env bash

set -e

configfile=$1
eval "$(egrep '^#|^[^ ]*=[^;&]*'  "$configfile")"
########
rm -rf ${output_folder}
mkdir ${output_folder}

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
java -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.PathSettingsCLI -temp_folder ${temp_folder} -log ${log_folder} &
wait

# create target-decoy fasta in the output directory
java -cp ${output_folder}/sg/SearchGUI-*.jar \
eu.isas.searchgui.cmd.FastaCLI \
-in ${fasta} -decoy &
wait

# Create search parameters
java -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -out ${output_folder}/identification.par -db ${fasta%.fasta}_concatenated_target_decoy.fasta -prec_tol ${prec_tol} -prec_ppm ${prec_ppm} -frag_tol ${frag_tol} -frag_ppm ${frag_ppm} -enzyme "${enzyme}" -fixed_mods "${fixed_mods}" -variable_mods "${variable_mods}" -min_charge ${min_charge} -max_charge ${max_charge} -mc ${mc} -fi ${fi} -ri ${ri} -psm_fdr ${psm_fdr} -peptide_fdr ${peptide_fdr} -protein_fdr ${protein_fdr} -myrimatch_min_pep_length ${myrimatch_min_pep_length} -myrimatch_max_pep_length ${myrimatch_max_pep_length} -msgf_instrument ${msgf_instrument} -msgf_min_pep_length ${msgf_min_pep_length} -msgf_max_pep_length ${msgf_max_pep_length} -tide_min_pep_length ${tide_min_pep_length} -tide_max_pep_length ${tide_max_pep_length} -import_peptide_length_min ${import_peptide_length_min} -import_peptide_length_max ${import_peptide_length_max} -annotation_level ${annotation_level} & 
wait

# database search function with three inputs:
function shakemgf {
    set -e
    spectrum_file=$1 #path to the spectrum
    experiment=$2    #name of the experiment
    replicate=$3     #replicate_number
    sample=$(basename $spectrum_file)
    
    echo "Started process for: "${sample}"..."
    sample_name=${sample%.mgf}
    sample_folder=${output_folder}"/"${sample_name}

    mkdir ${sample_folder} 
    mkdir ${sample_folder}/recalibrated
    mkdir ${sample_folder}/log
    mkdir ${sample_folder}/temp

    cp -R ${output_folder}/ps ${sample_folder}/ps
    cp -R ${output_folder}/sg ${sample_folder}/sg
    wait
	
    java -cp ${sample_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.PathSettingsCLI -temp_folder ${sample_folder}/temp -log ${sample_folder}/log  &
    wait

    echo "Started SearchCLI for: "${sample}"..." 
    java -cp ${sample_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.SearchCLI -spectrum_files ${spectrum_file} -output_folder ${sample_folder} -id_params ${output_folder}"/identification.par" -output_data "${output_data}" -xtandem "${xtandem}" -myrimatch "${myrimatch}" -ms_amanda "${ms_amanda}" -msgf "${msgf}" -omssa "${omssa}" -comet "${comet}" -tide "${tide}" -andromeda "${andromeda}" >> ${sample_folder}/${sample_name}.txt 2>&1 & 
    wait

    echo "Started PeptideShakerCLI for: "${sample}"..." 
    java -cp ${sample_folder}/ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI -experiment ${experiment} -sample ${sample} -replicate  ${replicate} -identification_files  ${sample_folder}"/searchgui_out.zip"  -spectrum_files  ${spectrum_file} -id_params  ${output_folder}"/identification.par" -out ${sample_folder}"/"${sample_name}".cpsx" >> ${sample_folder}/${sample_name}.txt 2>&1 & 
    wait

    echo "Started FollowUpCLI recalibration for: "${sample}"..." 
    java -cp ${sample_folder}/ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.FollowUpCLI -in ${sample_folder}"/"${sample_name}".cpsx" -recalibration_folder ${sample_folder}/recalibrated >> ${sample_folder}/${sample_name}.txt 2>&1 & 
    wait

    echo "Started recalibrated sample SearchCLI for: "${sample}"..." 
    java -cp ${sample_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.SearchCLI -spectrum_files ${sample_folder}/recalibrated/${sample_name}_recalibrated.mgf  -output_folder ${sample_folder}/recalibrated -id_params ${output_folder}"/identification.par" -output_data "${output_data}" -xtandem "${xtandem}" -myrimatch "${myrimatch}" -ms_amanda "${ms_amanda}" -msgf "${msgf}" -omssa "${omssa}" -comet "${comet}" -tide "${tide}" -andromeda "${andromeda}"  >> ${sample_folder}/${sample_name}.txt 2>&1 & 
    wait

    echo "Started recalibrated sample PeptideShakerCLI for: "${sample}"..."  
    java -cp ${sample_folder}/ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI -experiment ${experiment} -sample ${sample} -replicate ${replicate} -identification_files ${sample_folder}"/recalibrated/searchgui_out.zip" -spectrum_files ${sample_folder}/recalibrated/${sample_name}_recalibrated.mgf -id_params ${output_folder}"/identification.par" -out ${sample_folder}/recalibrated/${sample_name}".cpsx"  >> ${sample_folder}/${sample_name}.txt 2>&1 & 
    wait

    echo "Started recalibrated sample MzidCLI for: "${sample}"..." 
    java -cp ${sample_folder}/ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.MzidCLI -in ${sample_folder}/recalibrated/${sample_name}".cpsx" -output_file ${mzid_folder}"/"${sample_name}".mzid" -contact_first_name "${contact_first_name}" -contact_last_name "${contact_last_name}" -contact_email "${contact_email}" -contact_address "${contact_address}" -organization_name "${organization_name}" -organization_email "${organization_email}" -organization_address "${organization_address}" -contact_url "${contact_url}" -organization_url "${organization_url}" >> ${sample_folder}/${sample_name}.txt 2>&1 & 
    wait

    echo "Started recalibrated sample FollowUpCLI to export unvalidated spectra for: "${sample}"..."  
    java -cp ${sample_folder}/ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.FollowUpCLI -in ${sample_folder}/recalibrated/${sample_name}.cpsx -spectrum_folder ${unvalidated_folder} -psm_type ${psm_type} >> ${sample_folder}/${sample_name}.txt 2>&1 & 
    wait

    rm -rf ${sample_folder}/ps
    rm -rf ${sample_folder}/sg
}

count=0
for spectrum_file in ${spectrum_files}/*.mgf
do
  replicate=${count}
  shakemgf ${spectrum_file} ${experiment_name} ${replicate} & 
  let count+=1
done 
wait # wait for all the process to finish

echo 'MSnID FDR control' 
cd ${output_folder} && MSnIDshake.R -i mzIdentMLs/ -v ${MSnID_FDR_value} -l ${MSnID_FDR_level}

