#!/usr/bin/env bash

set -e
set -a
echo "MetaNovo version 9.0"
res1=$(date +%s.%N)
source /root/config.sh
source /root/bin/bash/compomics.sh

JVM_ARGS="-d64 -Xms${JVM_Xms} -Xmx${JVM_Xmx} -server"

mgf_file_count=$( find /root/mgf -name "*.mgf" | wc -l  )

 # Check that config does not exist or is unchanged
if [ ! -f /root/output/config.sh ] ; then
    cp /root/config.sh /root/output/config.sh
else
    cmp --silent /root/config.sh /root/output/config.sh && echo "'${CONFIG_FILE}' unchanged."|| { echo "'${CONFIG_FILE}' has changed, please delete '${OUTPUT_FOLDER}' or replace '${CONFIG_FILE}' with the contents of config.sh in "${OUTPUT_FOLDER}; exit 1; }
fi

if [ "$mn_search_database" -eq "1" ] ; then
    if [ ! -f /root/output/default_input.xml ] ; then
        cp /root/default_input.xml /root/output/default_input.xml
        cp /root/tandem-input-style.xsl /root/output/tandem-input-style.xsl
        echo "Please edit X!Tandem default_input.xml and tandem-input-style.xsl in output folder and restart the pipeline" && exit 1
    fi
fi

if [ ! -d /root/output/sg ] ; then
    cp -R SearchGUI* /root/output/sg
fi

if [ ! -d /root/output/utilities ] ; then
    cp -R utilities* /root/output/utilities
fi

output_folder=/root/output
input_fasta=/root/$(basename ${FASTA_FILE})

if [ ! -d ${output_folder}/denovo ] ; then
    mkdir ${output_folder}/denovo
    mkdir ${output_folder}/denovo/mgf 
    cp /root/mgf/*.mgf ${output_folder}/denovo/mgf
    mkdir ${output_folder}/denovo/runs
    mkdir ${output_folder}/denovo/temp
    mkdir ${output_folder}/denovo/log
fi

if [ ! -d /root/output/denovo/dg ] ; then
    cp -R DeNovoGUI* /root/output/denovo/dg
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

find ${output_folder}/denovo/mgf -name "*.mgf" \
    | parallel denovogui ${output_folder}/denovo {} ${output_folder}/identification.par $tagdb || exit 1

find ${output_folder}/denovo/mgf -name "*.gz"
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

find ${output_folder}/db -name "*.fasta" \
    | parallel -j${THREAD_LIMIT} "java -cp ${output_folder}/utilities/utilities-*.jar com.compomics.util.experiment.identification.protein_inference.executable.PeptideMapping -t {} ${output_folder}/tags/tags.txt {}.csv ${output_folder}/identification.par && bp_mapped_tags.py {}.csv {} ${tagdb} && gzip --best {}"

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
       mkdir $output_folder/mgf && cp /root/mgf/*.mgf ${output_folder}/mgf || { rm -rf $output_folder/mgf/ && echo "error copying mgf files"; exit 1; } 
    fi
    find ${output_folder}/mgf -name "*.mgf" \
        | parallel -j${THREAD_LIMIT} "rm -rf {}.*.xml && xtandem.R --mgf {} --fasta $input_fasta --input_xml $output_folder/default_input.xml --input_xsl $output_folder/tandem-input-style.xsl --output_xml {}.xml && gzip --best {}"
    
    find ${output_folder}/mgf -name "*.xml" \
        | parallel -j${THREAD_LIMIT} "java -Xms1024m -jar /mzidlib-*/mzidlib-*.jar Tandem2mzid {} {}.mzid -outputFragmentation false -idsStartAtZero false -decoyRegex :reversed -massSpecFileFormatID MS:1001062 -databaseFileFormatID MS:1001348 && gzip --best {}"
    if [ ! -d ${output_folder}/mgf/analysis ] ; then
        cd ${output_folder} && msgf_msnid.R -i mgf -v 1 -l "accession"  || rm -rf ${output_folder}/mgf/analysis
    fi
    if [ ! -f ${output_folder}/mgf/analysis/pept2lca.csv ] ; then 
        cd ${output_folder}/mgf/analysis && unipept pept2lca -i peptides_cleaned.txt -e -o pept2lca.csv || \
            rm -rf ${output_folder}/mgf/analysis/pept2lca.csv
    fi
    cd $output_folder/mgf/analysis && unipept.R
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

