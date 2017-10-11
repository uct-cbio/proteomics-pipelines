#!/usr/bin/env bash
set -e
sudo ./build.sh
txt_folder=/home/thys/bio/metanovo_stool_txt
#txt_folder=/home/thys/bio/stool/cns_reference
sudo docker run -it --rm -v ${txt_folder}:/root/txt cbio/mqmetaproteomics:latest

# Load data to elasticsearch
#cd ../csv2es && sudo ./run.sh 192.168.7.4 9230 "s507_s5527" "elastic" "changeme" ${outfolder}

