#!/usr/bin/env bash

set -e
JVM_ARGS="-d64 -Xms1024M -Xmx15360M -server"
configfile=$1
eval "$(egrep '^#|^[^ ]*=[^;&]*'  "$configfile")"
########
rm -rf ${output_folder}
mkdir ${output_folder}

temp_folder=${output_folder}"/temp"
log_folder=${output_folder}"/log"

mkdir ${temp_folder}
mkdir ${log_folder}

cp -R ${dg_folder} ${output_folder}/dg
wait


# Create search parameters

java $JVM_ARGS -cp ${output_folder}/dg/DeNovoGUI-*.jar com.compomics.denovogui.cmd.IdentificationParametersCLI -out ${output_folder}/identification.par -fixed_mods "${fixed_mods}" -variable_mods "${variable_mods}" -frag_tol ${frag_tol}   
wait

java $JVM_ARGS -cp ${output_folder}/dg/DeNovoGUI-*.jar com.compomics.denovogui.cmd.PathSettingsCLI -temp_folder ${temp_folder} -log ${log_folder}  
wait
java $JVM_ARGS -cp ${output_folder}/dg/DeNovoGUI-*.jar  com.compomics.denovogui.cmd.DeNovoCLI -spectrum_files ${spectrum_files} -output_folder ${output_folder} -id_params ${output_folder}/identification.par
wait
rm -rf ${output_folder}/dg

