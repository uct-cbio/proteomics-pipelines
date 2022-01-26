#/usr/bin/env bash

if [ $# -eq 0 ]                                                    
then                                                             
echo "Needs a config file as parameter"                         
exit 1                                                         
fi                                                                 
set -a

metanovo_version=v1.9.2

CONFIG_FILE=$(realpath $1)
echo $CONFIG_FILE
. ${CONFIG_FILE}
FASTA_BASE=$(basename $FASTA_FILE)
repo_path=$HOME/repos/proteomics-pipelines
sudo docker run -it --rm -v ${MGF_FOLDER}:/root/my_metanovo_project/mgf_files \
    -v ${FASTA_FILE}:/root/my_metanovo_project/${FASTA_BASE} \
    -v ${OUTPUT_FOLDER}:/root/my_metanovo_project/output \
    -v ${CONFIG_FILE}:/root/my_metanovo_project/config.sh  \
    -v ${repo_path}/bin:/home/bio/bin \
    -v ${repo_path}/lib:/home/bio/lib \
    thyscbio/metanovo:${metanovo_version} metanovo.sh /root/my_metanovo_project/config.sh 
