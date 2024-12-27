#!/usr/bin/env bash
set -e

sudo ./build.sh

outfolder=/data/bio/mtb/proteogenomics/multistrain1_test
interproscan_path="/data/bio/interproscan/5_31_70/interproscan-5.31-70.0"
sudo docker run -it --rm -v ${outfolder}:/root/output \
    -v $HOME/repos/proteomics-pipelines/bin/config/mq_proteogeomics_example_config.yml:/root/config.yaml \
    -v /data/bio/mtb/mq/multistrain1_test/txt:/root/txt \
    -v /data/bio/gbrowse/data/gff3/Mycobacterium_tuberculosis_h37rv.ASM19595v2.34.chromosome.Chromosome.gff3:/root/features.gff3 \
    -v /data/bio/mtb/genomes:/root/genomes \
    -v /data/bio/mtb/operons/3354.opr:/root/3354.txt \
    -v ${interproscan_path}:/interproscan \
    cbio/mqproteogenomics:latest

# Load data to elasticsearch
#cd ../csv2es && sudo ./run.sh 192.168.7.4 9210 "test" "elastic" "changeme" ${outfolder}
 

