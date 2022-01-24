#!/usr/bin/env bash

set -e
set -a


res1=$(date +%s.%N)
if [ $# -eq 0 ]; then
            CONFPATH="$( cd "$(dirname "$0")" && cd ../config >/dev/null 2>&1 ; pwd -P )"
            echo $CONFPATH
            cat ${CONFPATH}/metanovo_config.sh; exit 0
fi

config_file=$1
if [ ! -f "$config_file" ]; then
      echo $config_file does not exist; exit 1
fi

source ${config_file}

mgf_folder=${MGF_FOLDER}
if [ ! -d "$mgf_folder" ]; then
      echo $mgf_folder does not exist; exit 1
fi

fasta_file=${FASTA_FILE}
if [ ! -f "$fasta_file" ]; then
      echo $fasta_file does not exist; exit 1
fi

output_folder=${OUTPUT_FOLDER}/metanovo

export TMPDIR=${output_folder}/temp

if [ ! -d $TMPDIR ]; then
    mkdir -p $TMPDIR
fi

if [ ! -d "$output_folder" ]; then
      echo $output_folder does not exist; exit 1
fi

echo ${config_file}
echo ${mgf_folder}
echo ${fasta_file}
echo ${output_folder}

# set the TMPDIR to outputfolder...SQLITE crashes on ilifu..TO DO


source compomics.sh

JVM_ARGS="-Xms${JVM_Xms} -Xmx${JVM_Xmx} -server"
mgf_file_count=$( find ${mgf_folder} -name "*.mgf" | wc -l  )
if [ ! -d ${output_folder} ] ; then
    mkdir $output_folder
fi

CONFIG_FILE=${output_folder}/config.sh
 # Check that config does not exist or is unchanged
if [ ! -f ${output_folder}/config.sh ] ; then
    cp ${config_file} ${CONFIG_FILE}
else
    cmp --silent ${config_file} ${output_folder}/config.sh && echo "'${CONFIG_FILE}' unchanged."|| { echo "'${CONFIG_FILE}' has changed, please delete '${OUTPUT_FOLDER}' or replace '${CONFIG_FILE}' with the contents of config.sh in "${OUTPUT_FOLDER}; exit 1; }
fi

input_fasta=$FASTA_FILE
source `which env_parallel.bash`

if [ ! -d ${output_folder}/sg ] ; then
    cp -R ${SG_PATH} ${output_folder}/sg
fi

if [ ! -d ${output_folder}/utilities ] ; then
    cp -R ${CU_PATH} ${output_folder}/utilities
fi
if [ ! -d ${output_folder}/denovo ] ; then
    mkdir ${output_folder}/denovo
    mkdir ${output_folder}/denovo/mgf 
    cp ${mgf_folder}/*.mgf ${output_folder}/denovo/mgf
    mkdir ${output_folder}/denovo/runs
    mkdir ${output_folder}/denovo/temp
    mkdir ${output_folder}/denovo/log
fi
if [ ! -d ${output_folder}/denovo/dg ] ; then
    cp -R ${DG_PATH} ${output_folder}/denovo/dg
fi
if [ ! -f ${output_folder}/identification.par ] ; then
    parameters
fi
denovo_file_count=$( find ${output_folder}/denovo/mgf -name "*.mgf*" -not -name "*.cui" | wc -l  )
if [ ! ${denovo_file_count} -eq ${mgf_file_count} ]; then
    echo "Not all mgfs have been copied" && exit 1
fi
tagdb=$output_folder/sqlite3.db
#export -f denovogui

# Denovo sequencing
cmd="denovogui ${output_folder}/denovo {} ${output_folder}/identification.par $tagdb || exit 1 "
if [ -z "${PBS_NODEFILE}" ] ; then
  find ${output_folder}/denovo/mgf -name "*.mgf" | parallel -q -j${THREAD_LIMIT} ${cmd}
else
  find ${output_folder}/denovo/mgf -name "*.mgf" | env_parallel -q -j${THREAD_LIMIT} --sshloginfile ${PBS_NODEFILE} ${cmd}
fi

if [ ! -d $output_folder/tags ] ; then
    mkdir $output_folder/tags
fi

if [ ! -f $output_folder/tags/tags.txt ] ; then
    bp_export_tags.py ${output_folder}/tags/tags.txt $tagdb || exit 1
fi
if [ ! -d ${output_folder}/db ] ; then
    echo "PREPARING FASTA FOR SEARCHES"
    mkdir ${output_folder}/db && \
        bp_fasta_prepare.py ${input_fasta} ${CHUNKSIZE}  ${output_folder}/db/ \
        || { rm -rf ${output_folder}/db ; echo "Error in database splitting in analysis section";     exit 1; }
fi

cmd="java ${JVM_ARGS} -cp ${output_folder}/utilities/utilities-*.jar com.compomics.util.experiment.identification.protein_inference.executable.PeptideMapping -t {} ${output_folder}/tags/tags.txt {}.csv ${output_folder}/identification.par && gzip --best {}"

#cmd="java ${JVM_ARGS} -cp ${output_folder}/utilities/utilities-*.jar com.compomics.cli.peptide_mapper.PeptideMapperCLI -t {} ${output_folder}/tags/tags.txt {}.csv ${output_folder}/identification.par && gzip --best {}"

if [ -z "${PBS_NODEFILE}" ]; then
    find ${output_folder}/db -name "*.fasta" \
    | parallel -j${THREAD_LIMIT} ${cmd}
else
    find ${output_folder}/db -name "*.fasta" \
    | env_parallel -j${THREAD_LIMIT} --sshloginfile ${PBS_NODEFILE} ${cmd}
fi

cmd="bp_mapped_tags.py {} ${tagdb} && gzip --best {}"
if [ -z "${PBS_NODEFILE}" ]; then
    find ${output_folder}/db -name "*.csv" \
        | parallel -j${THREAD_LIMIT} ${cmd} 
else
    find ${output_folder}/db -name "*.csv" \
        | env_parallel -j${THREAD_LIMIT} --sshloginfile ${PBS_NODEFILE} ${cmd}
fi

if [ ! -f ${output_folder}/metanovo.fasta ] ; then
    bp_export_proteins.py ${tagdb} ${output_folder}
fi

###################
# time the script #
###################
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)
printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds

