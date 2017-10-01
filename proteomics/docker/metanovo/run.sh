#/usr/bin/env bash

if [ $# -eq 0 ]                                                    
then                                                             
echo "Needs a config file as parameter"                         
exit 1                                                         
fi                                                                 
set -a

CONFIG_FILE=$(realpath $1)

source ${CONFIG_FILE}

FASTA_BASE=$(basename $FASTA_FILE)

sudo docker run -it --rm -v ${MGF_FOLDER}:/root/mgf \
    -v ${FASTA_FILE}:/root/${FASTA_BASE} \
    -v ${OUTPUT_FOLDER}:/root/output \
    -v ${CONFIG_FILE}:/root/config.sh  \
    -e CONFIG_FILE=${CONFIG_FILE} \
    cbio/metanovo:latest metanovo.sh /root/mgf /root/${FASTA_BASE} /root/output /root/config.sh 
