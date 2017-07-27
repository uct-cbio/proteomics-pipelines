#!/usr/bin/env bash

set -e
output_folder=$1
spectrum_file=$2 #path to the spectrum

JVM_ARGS="-d64 -Xms1024M -Xmx15360M -server"
configfile=$output_folder/config.cfg
eval "$(egrep '^#|^[^ ]*=[^;&]*'  "$configfile")"
sample=$(basename $spectrum_file)

echo "Started process for: "${sample}"..."
sample_name=${sample%.mgf}
sample_folder=${output_folder}"/"${sample_name}
replicate=0

mzid_folder=${output_folder}"/mzIdentMLs"


rm -rf $sample_folder
mkdir ${sample_folder} 
mkdir ${sample_folder}/log
mkdir ${sample_folder}/temp

cp -R ${output_folder}/ps ${sample_folder}/ps
cp -R ${output_folder}/sg ${sample_folder}/sg

ps=${sample_folder}/ps
sg=${sample_folder}/sg
log=${sample_folder}/log/${sample_name}.mgf.out

wait

java $JVM_ARGS -cp ${sample_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.PathSettingsCLI -temp_folder ${sample_folder}/temp -log ${sample_folder}/log  ; wait

function searchgui {	
    java $JVM_ARGS -cp $sg/SearchGUI-*.jar eu.isas.searchgui.cmd.SearchCLI -spectrum_files ${spectrum_file} -output_folder ${sample_folder} -id_params ${output_folder}"/identification.par" -output_data "${output_data}" -xtandem "${xtandem}" -myrimatch "${myrimatch}" -ms_amanda "${ms_amanda}" -msgf "${msgf}" -omssa "${omssa}" -comet "${comet}" -tide "${tide}" -andromeda "${andromeda}" >> $log 2>&1 ; wait 
}

function peptideshaker {
    java $JVM_ARGS -cp $ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI -experiment ${experiment} -sample $sample -replicate  ${replicate} -identification_files  ${sample_folder}"/searchgui_out.zip"  -spectrum_files  ${spectrum_file} -id_params  ${output_folder}"/identification.par" -out ${sample_folder}"/"${sample}".cpsx" >> $log ; wait 
}

function recalibrate {
    java $JVM_ARGS -cp $ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.FollowUpCLI -in ${sample_folder}"/"$sample".cpsx" -recalibration_folder ${sample_folder}/recalibrated >> $log 2>&1 ; wait 
}

function mzidcli {
    java $JVM_ARGS -cp $ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.MzidCLI -in ${sample_folder}/$sample".cpsx" -output_file ${mzid_folder}"/"${sample}".mzid" -contact_first_name "${contact_first_name}" -contact_last_name "${contact_last_name}" -contact_email "${contact_email}" -contact_address "${contact_address}" -organization_name "${organization_name}" -organization_email "${organization_email}" -organization_address "${organization_address}" -contact_url "${contact_url}" -organization_url "${organization_url}" >> $log 2>&1 ; wait 
}

function unvalidated {
    java $JVM_ARGS -cp $ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.FollowUpCLI -in ${sample_folder}/${sample}.cpsx -spectrum_folder ${unvalidated_folder} -psm_type ${psm_type}>> $log 2>&1 ; wait 
}

searchgui && peptideshaker

if [ $recalibrate -eq 1 ]; then   
    mkdir ${sample_folder}/recalibrated && recalibrate 
    spectrum_file=${sample_folder}/recalibrated/${sample_name}_recalibrated.mgf
    sample=$(basename $spectrum_file)
    log=${sample_folder}/log/${sample_name}_recalibrated.mgf.out
    sample_folder=${sample_folder}/recalibrated 
    searchgui && peptideshaker 
fi

mzidcli && unvalidated

rm -rf ${sample_folder}/ps
rm -rf ${sample_folder}/sg

