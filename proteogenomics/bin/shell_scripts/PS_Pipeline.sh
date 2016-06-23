#!/usr/bin/env bash

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


fasta=${fasta_folder}"/"$(basename ${target_fasta})
unvalidated_folder=${output_folder}"/unvalidated"

mkdir ${unvalidated_folder}

temp_folder=${output_folder}"/temp"
mkdir ${temp_folder}

cp -R ${ps_folder} ${temp_folder}/ps
cp -R ${sg_folder} ${temp_folder}/sg

# create target-decoy fasta in the output directory
#cd ${HOME}/SearchGUI-2.9.0
java -cp ${temp_folder}/sg/SearchGUI-*.jar \
eu.isas.searchgui.cmd.FastaCLI \
-in ${fasta} -decoy


# create search parameters
java -cp ${temp_folder}/ps/PeptideShaker-*.jar \
eu.isas.peptideshaker.cmd.IdentificationParametersCLI \
-out                 "${output_folder}""/identification.par" \
-db                  "${fasta%.fasta}""_concatenated_target_decoy.fasta" \
-prec_tol            "${prec_tol}" \
-prec_ppm            "${prec_ppm}" \
-frag_tol            "${frag_tol}" \
-frag_ppm            "${frag_ppm}" \
-enzyme              "${enzyme}"   \
-fixed_mods          "${fixed_mods}" \
-variable_mods       "${variable_mods}" \
-min_charge          "${min_charge}" \
-max_charge          "${max_charge}" \
-fi                  "${fi}" \
-ri                  "${ri}"\
-mc		     "${mc}"	

# database search function with three inputs:
function shakemgf {
    spectrum_file=$1 #path to the spectrum
    experiment=$2    #name of the experiment
    replicate=$3     #replicate_number

    sample=$(basename $spectrum_file)
    sample_name=${sample%.mgf}
    sample_folder=${output_folder}"/"${sample_name}
    mkdir ${sample_folder} 

    java -cp ${temp_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.SearchCLI \
    -spectrum_files      ${spectrum_file} \
    -output_folder       ${sample_folder}  \
    -id_params           ${output_folder}"/identification.par" \
    -output_data         "${output_data}" \
    -xtandem             "${xtandem}" \
    -myrimatch           "${myrimatch}" \
    -ms_amanda           "${ms_amanda}" \
    -msgf                "${msgf}" \
    -omssa               "${omssa}" \
    -comet               "${comet}" \
    -tide                "${tide}" \
    -andromeda           "${andromeda}"

    java -cp ${temp_folder}/ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI \
    -experiment            ${experiment} \
    -sample                ${sample} \
    -replicate             ${replicate} \
    -identification_files  ${sample_folder}"/searchgui_out.zip"  \
    -spectrum_files        ${spectrum_file} \
    -id_params             ${output_folder}"/identification.par" \
    -out                   ${sample_folder}"/"${sample_name}".cpsx"

    java -cp ${temp_folder}/ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.MzidCLI \
    -in                    ${sample_folder}"/"${sample_name}".cpsx" \
    -output_file           ${mzid_folder}"/"${sample_name}".mzid" \
    -contact_first_name    "${contact_first_name}" \
    -contact_last_name     "${contact_last_name}" \
    -contact_email         "${contact_email}" \
    -contact_address       "${contact_address}" \
    -organization_name     "${organization_name}" \
    -organization_email    "${organization_email}" \
    -organization_address  "${organization_address}" \
    -contact_url           "${contact_url}" \
    -organization_url      "${organization_url}"

    java -cp ${temp_foler}/ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.FollowUpCLI \
    -in                     ${sample_folder}"/"${sample_name}".cpsx" \
    -spectrum_folder        ${unvalidated_folder} \
    -psm_type               "${psm_type}"
}

count=0
for spectrum_file in ${spectrum_files}/*.mgf
do
  echo ${spectrum_file}
  replicate=${count}
  shakemgf ${spectrum_file} ${experiment_name} ${replicate}
  let count+=1
done

cd ${output_folder} 

MSnIDshake.R -i mzIdentMLs/ --psm_fdr ${psm_FDR} --peptide_fdr {$pep_FDR} --protein_fdr {$prot_FDR}

