#!/usr/bin/env bash
set -e

function pathsettingscli {
java $JVM_ARGS -cp ${sample_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.PathSettingsCLI -temp_folder ${sample_folder}/temp -log ${sample_folder}/log  ; wait
}

function searchgui {	
    java $JVM_ARGS -cp $sg/SearchGUI-*.jar eu.isas.searchgui.cmd.SearchCLI -mgf_splitting "${mgf_splitting}" -threads ${threads} -spectrum_files ${spectrum_file} -output_folder ${sample_folder} -id_params ${output_folder}"/identification.par" -output_data "${output_data}" -xtandem "${xtandem}" -myrimatch "${myrimatch}" -ms_amanda "${ms_amanda}" -msgf "${msgf}" -omssa "${omssa}" -comet "${comet}" -tide "${tide}" -andromeda "${andromeda}" >> $log 2>&1 ; wait 
}

function peptideshaker {
    java $JVM_ARGS -cp $ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI -experiment ${experiment} -threads ${threads} -sample $sample -replicate  ${replicate} -identification_files  ${sample_folder}"/searchgui_out.zip"  -spectrum_files  ${spectrum_file} -id_params  ${output_folder}"/identification.par" -out ${sample_folder}"/"${sample}".cpsx" >> $log ; wait 
}

function recalibrate {
    java $JVM_ARGS -cp $ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.FollowUpCLI -in ${sample_folder}"/"$sample".cpsx" -recalibration_folder ${sample_folder}/recalibrated >> $log 2>&1 ; wait 
}

function mzidcli {
    java $JVM_ARGS -cp $ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.MzidCLI -in ${sample_folder}/$sample".cpsx" -output_file ${sample_folder}"/"${sample}".mzid" -contact_first_name "${contact_first_name}" -contact_last_name "${contact_last_name}" -contact_email "${contact_email}" -contact_address "${contact_address}" -organization_name "${organization_name}" -organization_email "${organization_email}" -organization_address "${organization_address}" -contact_url "${contact_url}" -organization_url "${organization_url}" >> $log 2>&1 && mv ${sample_folder}"/"${sample}".mzid" ${mzid_folder} && gzip --best $original_spectrum_file  ; wait 
}

function unvalidated {
    java $JVM_ARGS -cp $ps/PeptideShaker-*.jar eu.isas.peptideshaker.cmd.FollowUpCLI -in ${sample_folder}/${sample}.cpsx -spectrum_folder ${unvalidated_folder} -psm_type ${psm_type}>> $log 2>&1 ; wait 
}

function ps_prepare {
    set -e
    mzid_folder=${output_folder}"/searchgui/mzIdentMLs"
    mkdir ${mzid_folder} || { rm -rf $output_folder/searchgui ; exit 1 ; }

    fasta_folder=${output_folder}"/searchgui/fasta"
    mkdir ${fasta_folder} || rm -rf $output_folder/searchgui

    mgf_folder=${output_folder}"/searchgui/mgf"
    mkdir ${mgf_folder} || rm -rf $output_folder/searchgui

    cp ${target_fasta} ${fasta_folder}
    
    
    cp ${spectrum_files}/*.mgf $mgf_folder || rm -rf $output_folder/searchgui

    target=$(basename $target_fasta)
    
    target_fasta=${fasta_folder}/${target}

    target=${target%.fasta}

    unvalidated_folder=${output_folder}"/searchgui/unvalidated"

    mkdir ${unvalidated_folder}

    temp_folder=${output_folder}"/searchgui/temp"
    log_folder=${output_folder}"/searchgui/log"

    mkdir ${temp_folder} || { rm -rf $output_folder/searchgui ; exit 1 ; }
    mkdir ${log_folder} || { rm -rf $output_folder/searchgui ; exit 1 ; }
    mkdir ${output_folder}/searchgui/runs || { rm -rf $output_folder/searchgui ; exit 1 ; }
    wait
    # Set paths for logs and temporary files
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.PathSettingsCLI -temp_folder ${temp_folder} -log ${log_folder} 
    wait
    # create target-decoy fasta in the output directory
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar \
    eu.isas.searchgui.cmd.FastaCLI \
    -in ${target_fasta} -decoy 
    wait
    # Create search parameters
    # Set FASTA file
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -out ${output_folder}/identification.par -db ${target_fasta%.fasta}_concatenated_target_decoy.fasta 
    wait

    }


function search {
    output_folder=$1
    spectrum_file=$2 
    sample=$(basename $spectrum_file)
    original_spectrum_file=$spectrum_file
    echo "Started process for: "${sample}"..."
    sample_name=${sample%.mgf}
    sample_folder=${output_folder}"/searchgui/runs/"${sample_name}
    replicate=0
    mzid_folder=${output_folder}"/searchgui/mzIdentMLs"
    unvalidated_folder=${output_folder}"/searchgui/unvalidated"
    rm -rf $sample_folder
    mkdir ${sample_folder} 
    mkdir ${sample_folder}/log
    mkdir ${sample_folder}/temp
    cp -R ${output_folder}/ps ${sample_folder}/ps
    cp -R ${output_folder}/sg ${sample_folder}/sg
    ps=${sample_folder}/ps
    sg=${sample_folder}/sg
    log=${sample_folder}/log/${sample_name}.mgf.out
    pathsettingscli
    searchgui && peptideshaker
    unvalidated && mzidcli
    rm -rf ${sample_folder}/ps
    rm -rf ${sample_folder}/sg
    wait
}


function parameters {
    # Create search parameters
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -out ${output_folder}/identification.par  \
        \
        -prec_tol ${prec_tol} \
        -prec_ppm ${prec_ppm} \
        -frag_tol ${frag_tol} \
        -frag_ppm ${frag_ppm} \
        -enzyme "${enzyme}" \
        -specificity "${specificity}" \
        -digestion "${digestion}" \
        -fixed_mods "${fixed_mods}" \
        -variable_mods "${variable_mods}" \
        -min_charge ${min_charge} \
        -max_charge ${max_charge} \
        -min_isotope ${min_isotope} \
        -max_isotope ${max_isotope} \
        -mc ${mc} \
        -fi ${fi} \
        -ri ${ri} \
        \
        -annotation_level ${annotation_level}  \
        \
        -sequence_index_type ${sequence_index_type} \
        -sequence_matching_type ${sequence_matching_type} \
        -sequence_matching_x ${sequence_matching_x} \
        \
        -import_peptide_length_min ${import_peptide_length_min} -import_peptide_length_max ${import_peptide_length_max} -import_precurosor_mz_ppm ${import_precursor_mz_ppm} -exclude_unknown_ptms ${exclude_unknown_ptms} \
       -ptm_score ${ptm_score} -score_neutral_losses ${score_neutral_losses} -ptm_sequence_matching_type ${ptm_sequence_matching_type} -ptm_alignment ${ptm_alignment} \
        -useGeneMapping ${useGeneMapping} \
        -updateGeneMapping ${updateGeneMapping} \
        \
        -simplify_groups ${simplify_groups} \
        -simplify_score ${simplify_score} \
        -simplify_enzymaticity ${simplify_enzymaticity} \
        -simplify_evidence ${simplify_evidence} \
        -simplify_uncharacterized ${simplify_uncharacterized} \
        -psm_fdr ${psm_fdr} \
        -peptide_fdr ${peptide_fdr} \
        -protein_fdr ${protein_fdr} \
        -group_psms ${group_psms} \
        -group_peptides ${group_peptides} \
        -merge_subgroups ${merge_subgroups} \
        \
        -protein_fraction_mw_confidence ${protein_fraction_mw_confidence} \
        \
        -pepnovo_hitlist_length ${pepnovo_hitlist_length} \
        -pepnovo_estimate_charge ${pepnovo_estimate_charge} \
        -pepnovo_correct_prec_mass ${pepnovo_correct_prec_mass} \
        -pepnovo_discard_spectra ${pepnovo_discard_spectra} \
        -pepnovo_fragmentation_model ${pepnovo_fragmentation_model} \
        -pepnovo_generate_blast ${pepnovo_generate_blast} \
        \
        -directag_tic_cutoff ${directag_tic_cutoff} \
        -directag_max_peak_count ${directag_max_peak_count} \
        -directag_intensity_classes ${directag_intensity_classes} \
        -directag_adjust_precursor ${directag_adjust_precursor} \
        -directag_min_adjustment ${directag_min_adjustment} \
        -directag_max_adjustment ${directag_max_adjustment} \
        -directag_adjustment_step ${directag_adjustment_step} -directag_charge_states ${directag_charge_states} -directag_ms_charge_state ${directag_ms_charge_state} -directag_duplicate_spectra ${directag_duplicate_spectra} -directag_deisotoping ${directag_deisotoping} -directag_isotope_tolerance ${directag_isotope_tolerance} -directag_complement_tolerance ${directag_complement_tolerance} -directag_tag_length ${directag_tag_length} -directag_max_var_mods ${directag_max_var_mods} -directag_max_tag_count ${directag_max_tag_count} -directag_intensity_weight ${directag_intensity_weight} -directag_fidelity_weight ${directag_fidelity_weight} -directag_complement_weight ${directag_complement_weight} \
        -novor_fragmentation ${novor_fragmentation} -novor_mass_analyzer ${novor_mass_analyzer}
    wait

    }




function denovogui {
    set -e
    output_folder="$1"
    spectrum_file="$2" 
    id_params="$3"
    tagdb="$4"
    sample="$(basename "$spectrum_file")"
    original_spectrum_file="$spectrum_file"

    echo "Started process for: "${sample}"..."
    sample_name="${sample%.mgf}"
    sample_folder=${output_folder}"/runs/""${sample_name}"
    replicate=0

    rm -rf "$sample_folder"
    mkdir "${sample_folder}" 
    cp "${spectrum_file}" "${sample_folder}"

    mkdir "${sample_folder}/log"
    mkdir "${sample_folder}/temp"

    cp -R ${output_folder}/dg "${sample_folder}/dg"

    java -cp "${sample_folder}"/dg/DeNovoGUI-*.jar com.compomics.denovogui.cmd.PathSettingsCLI -temp_folder "${sample_folder}/temp" \
        -log "${sample_folder}/log"
   

    java -cp "${sample_folder}"/dg/DeNovoGUI-*.jar com.compomics.denovogui.cmd.DeNovoCLI \
        -pepnovo ${dg_pepnovo} \
        -directag ${dg_directag} \
        -pnovo ${dg_pnovo} \
        -novor ${dg_novor} \
        -spectrum_files "${sample_folder}/${sample}" \
        -output_folder "${sample_folder}" \
        -id_params ${id_params} >> "${sample_folder}/log/${sample_name}.log" 2>&1 \
        && bp_parse_tags.py "${sample_folder}" "${sample}" ${tagdb} \
        && rm -rf "${sample_folder}/dg" \
        && tar -czf "${spectrum_file}.tar.gz" "${spectrum_file}"  \
        && rm -rf "${spectrum_file}" \
        && rm -rf "${sample_folder}/*.mgf" || exit 1
    wait
}
