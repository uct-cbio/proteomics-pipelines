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
    mzid_folder=${output_folder}"/mzIdentMLs"
    mkdir ${mzid_folder}

    fasta_folder=${output_folder}"/fasta"
    mkdir ${fasta_folder}

    mgf_folder=${output_folder}"/mgf"
    mkdir ${mgf_folder}

    cp ${target_fasta} ${fasta_folder}
    cp ${contaminant_fasta} ${fasta_folder}
    cp ${spectrum_files}/*.mgf $mgf_folder
    
    contaminant=$(basename $contaminant_fasta)
    target=$(basename $target_fasta)

    contaminant=${contaminant%.fasta}
    target=${target%.fasta}

    fasta=${fasta_folder}/${target}_${contaminant}.fasta
    cat ${target_fasta} <(annotateContaminantUniprotFasta.py ${contaminant_fasta}) > ${fasta}
    #fasta=${fasta_folder}/${target}.fasta 

    unvalidated_folder=${output_folder}"/unvalidated"

    mkdir ${unvalidated_folder}

    temp_folder=${output_folder}"/temp"
    log_folder=${output_folder}"/log"
    mkdir ${temp_folder}
    mkdir ${log_folder}
    mkdir ${output_folder}/runs
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
    # Set FASTA file
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -out ${output_folder}/identification.par -db ${fasta%.fasta}_concatenated_target_decoy.fasta 
    wait
    # Set Tolerances
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -prec_tol ${prec_tol} -prec_ppm ${prec_ppm} -frag_tol ${frag_tol} -frag_ppm ${frag_ppm} 
    wait
    # Set Enzyme
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -enzyme "${enzyme}" 
    wait
    # Set Modifications
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -fixed_mods "${fixed_mods}" -variable_mods "${variable_mods}"
    wait
    # Set Charges
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -min_charge ${min_charge} -max_charge ${max_charge} 
    wait
    # Set MC
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -mc ${mc} 
    wait
    # Set ion type
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -fi ${fi} -ri ${ri} 
    wait
    # Set FDR
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -psm_fdr ${psm_fdr} -peptide_fdr ${peptide_fdr} -protein_fdr ${protein_fdr} 
    wait
    # Myrimatch settings
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -myrimatch_min_pep_length ${myrimatch_min_pep_length} -myrimatch_max_pep_length ${myrimatch_max_pep_length} 
    wait
    # MSGF+ settings
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -msgf_instrument ${msgf_instrument} -msgf_min_pep_length ${msgf_min_pep_length} -msgf_max_pep_length ${msgf_max_pep_length} 
    wait
    
    # Tide settings
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -tide_min_pep_length ${tide_min_pep_length} -tide_max_pep_length ${tide_max_pep_length} 
    wait
    
    # Import settings
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -import_peptide_length_min ${import_peptide_length_min} -import_peptide_length_max ${import_peptide_length_max} 
    wait
    
    # Gene Annotation
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -useGeneMapping ${useGeneMapping} -updateGeneMapping ${updateGeneMapping}
    wait
    
    # Set annotation level
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -annotation_level ${annotation_level}
    wait

    # PTM localization
    java $JVM_ARGS -cp ${output_folder}/sg/SearchGUI-*.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -id_params ${output_folder}/identification.par -out ${output_folder}/identification.par -ptm_score ${ptm_score} -score_neutral_losses ${score_neutral_losses} -ptm_sequence_matching_type ${ptm_sequence_matching_type} -ptm_alignment ${ptm_alignment}
    wait

    }

function search {
    output_folder=$1
    spectrum_file=$2 
    sample=$(basename $spectrum_file)
    original_spectrum_file=$spectrum_file

    echo "Started process for: "${sample}"..."
    sample_name=${sample%.mgf}
    sample_folder=${output_folder}"/runs/"${sample_name}
    replicate=0

    mzid_folder=${output_folder}"/mzIdentMLs"
    unvalidated_folder=${output_folder}"/unvalidated"

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

    if [ $recalibrate -eq 1 ]; then
       mkdir ${sample_folder}/recalibrated && recalibrate
       spectrum_file=${sample_folder}/recalibrated/${sample_name}_recalibrated.mgf
       sample=$(basename $spectrum_file)
       log=${sample_folder}/log/${sample_name}_recalibrated.mgf.out
       sample_folder=${sample_folder}/recalibrated
       searchgui && peptideshaker
    fi
  
    unvalidated && mzidcli
  
    rm -rf ${sample_folder}/ps
    rm -rf ${sample_folder}/sg
  
    wait
}

