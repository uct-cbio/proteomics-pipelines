#!/usr/bin/env bash

set -e
set -a
echo "MetaNovo version 9.0"
res1=$(date +%s.%N)


mgf_folder=$1
fasta_file=$2
output_folder=$3/metanovo
config_file=$4

source ${config_file}

source compomics.sh

JVM_ARGS="-d64 -Xms${JVM_Xms} -Xmx${JVM_Xmx} -server"
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

if [ "$mn_search_database" -eq "1" ] ; then
    if [ ! -f ${output_folder}/default_input.xml ] ; then
        cp ${TANDEM_DEFAULT_INPUT_PATH} ${output_folder}/default_input.xml
        cp ${TANDEM_INPUT_STYLE_PATH} ${output_folder}/tandem-input-style.xsl
        echo "Please edit X!Tandem default_input.xml and tandem-input-style.xsl in OUTPUT_FOLDER/metanovo and restart the pipeline" && exit 1
    fi
fi

if [ ! -d ${output_folder}/sg ] ; then
    cp -R ${SG_PATH} ${output_folder}/sg
fi

if [ ! -d ${output_folder}/utilities ] ; then
    cp -R ${CU_PATH} ${output_folder}/utilities
fi


input_fasta=$2

if [ ! -d ${output_folder}/denovo ] ; then
    mkdir ${output_folder}/denovo
    mkdir ${output_folder}/denovo/mgf 
    cp ${mgf_folder}/*.mgf ${output_folder}/denovo/mgf
    mkdir ${output_folder}/denovo/runs
    mkdir ${output_folder}/denovo/temp
    mkdir ${output_folder}/denovo/log
fi

if [ ! -d ${output_fodler}/denovo/dg ] ; then
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
  find ${output_folder}/denovo/mgf -name "*.mgf" | parallel -j${THREAD_LIMIT} ${cmd}
else
  find ${output_folder}/denovo/mgf -name "*.mgf" | parallel -j${THREAD_LIMIT} --sshloginfile ${PBS_NODEFILE} ${cmd}
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
if [ -z "${PBS_NODEFILE}" ]; then
    find ${output_folder}/db -name "*.fasta" \
    | parallel -j${THREAD_LIMIT} ${cmd}
else
    find ${output_folder}/db -name "*.fasta" \
    | parallel -j${THREAD_LIMIT} --sshloginfile ${PBS_NODEFILE} ${cmd}
fi

cmd="bp_mapped_tags.py {} ${tagdb} && gzip --best {}"
if [ -z "${PBS_NODEFILE}" ]; then
    find ${output_folder}/db -name "*.csv" \
        | parallel -j${THREAD_LIMIT} ${cmd} 
else
    find ${output_folder}/db -name "*.csv" \
        | parallel -j${THREAD_LIMIT} --sshloginfile ${PBS_NODEFILE} ${cmd}
fi

if [ ! -f ${output_folder}/metanovo.fasta ] ; then
    bp_export_proteins.py ${tagdb} ${output_folder}
fi

function msgfplus {
    echo "Starting msgfplus"
    java ${JVM_ARGS} -jar /msgfplus/MSGFPlus.jar -s $1 -d $2 -o $3 -t $msgfplus_t -ti $msgfplus_ti -tda $msgfplus_tda -m $msgfplus_m -inst $msgfplus_inst -e $msgfplus_e -protocol $msgfplus_protocol -ntt $msgfplus_ntt -mod ${output_folder}/MSGFPlus_Mods.txt -minLength $msgfplus_minLength -maxLength $msgfplus_maxLength -minCharge $msgfplus_minCharge -maxCharge  $msgfplus_maxCharge -n $msgfplus_n -addFeatures $msgfplus_addFeatures -ccm $msgfplus_ccm
}

if [ "$mn_filter_database" -eq "1" ] ; then
    input_fasta=$output_folder/metanovo.fasta
fi

if [ "$mn_search_database" -eq "1" ] ; then
    if [ ! -d $output_folder/mgf ] ; then
       mkdir $output_folder/mgf && cp ${mgf_folder}/*.mgf ${output_folder}/mgf || { rm -rf $output_folder/mgf/ && echo "error copying mgf files"; exit 1; } 
    fi
    
    # Xtandem
    cmd="rm -rf {}.*.xml && xtandem.R --mgf {} --fasta $input_fasta --input_xml $output_folder/default_input.xml --input_xsl $output_folder/tandem-input-style.xsl --output_xml {}.xml && gzip --best {}" 
    if [ -z "${PBS_NODEFILE}" ]; then
        find ${output_folder}/mgf -name "*.mgf" \
            | parallel -j${THREAD_LIMIT} ${cmd}
    else
        find ${output_folder}/mgf -name "*.mgf" \
            | parallel -j${THREAD_LIMIT} --sshloginfile ${PBS_NODEFILE} ${cmd}
    fi
    
    # MZID
    cmd="java -Xms1024m -jar ${MZIDLIB_PATH}/mzidlib-*.jar Tandem2mzid {} {}.mzid -outputFragmentation false -idsStartAtZero false -decoyRegex :reversed -massSpecFileFormatID MS:1001062 -databaseFileFormatID MS:1001348 && gzip --best {}"
    if [ -z "${PBS_NODEFILE}" ]; then
        find ${output_folder}/mgf -name "*.xml" \
            | parallel -j${THREAD_LIMIT} ${cmd} 
    else
        find ${output_folder}/mgf -name "*.xml" \
            | parallel -j${THREAD_LIMIT} --shhloginfile ${PBS_NODEFILE} ${cmd}
    fi


    if [ ! -d ${output_folder}/mgf/analysis ] ; then
        cd ${output_folder} && msgf_msnid.R -i mgf -v ${mn_search_fdr_value} -l ${mn_search_fdr_level}  || rm -rf ${output_folder}/mgf/analysis
    fi

    if [ ! -f ${output_folder}/mgf/analysis/pept2lca.csv ] ; then 
        cd ${output_folder}/mgf/analysis && unipept pept2lca -i peptides_cleaned.txt -e -o pept2lca.csv || \
            rm -rf ${output_folder}/mgf/analysis/pept2lca.csv
    fi
    if [ ! -f ${output_folder}/mgf/analysis/accession_sc_pept2lca.txt  ] ; then
        cd $output_folder/mgf/analysis && unipept.R
    fi
    #if [ ! -d $output_folder/lfq ] ; then
    #   mkdir $output_folder/lfq && cp /root/mgf/*.mgf ${output_folder}/lfq || { rm -rf $output_folder/lfq/ && echo "error copying mgf files"; exit 1; } 
    #fi
    #cd ${output_folder} && find lfq -name "*.mgf" \
    #    | parallel -j${THREAD_LIMIT} "msnbase.R --mgf {} --psms ${output_folder}/mgf/analysis/tables/psms.txt"
fi


#if [ ! -f ${output_folder}/tags/trie.p ] ; then
#   create_trie.py $tagdb tags sequence /${output_folder}/tags/trie.p || { rm -rf ${output_folder}/tags/trie.p ; exit 1 ; }
#fi

#if [ ! -f ${output_folder}/tags/search.txt ] ; then
#    bp_triematching.py ${input_fasta} ${output_folder}/tags/trie.p ${tagdb} 
#fi


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

