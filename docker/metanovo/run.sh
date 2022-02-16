#/usr/bin/env bash

if [ $# -eq 0 ]                                                    
then                                                             
echo "Needs a config file as parameter"                         
exit 1                                                         
fi                                                                 
set -a

metanovo_version=v1.9.4

CONFIG_FILE=$(realpath $1)
echo $CONFIG_FILE
. ${CONFIG_FILE}
FASTA_BASE=$(basename $FASTA_FILE)
sudo docker run -it --rm -v ${MGF_FOLDER}:/root/my_metanovo_project/mgf_files \
    -v ${FASTA_FILE}:/root/my_metanovo_project/${FASTA_BASE} \
    -v ${OUTPUT_FOLDER}:/root/my_metanovo_project/output \
    -v ${CONFIG_FILE}:/root/my_metanovo_project/config.sh  \
    thyscbio/metanovo:${metanovo_version} metanovo.sh /root/my_metanovo_project/config.sh 
