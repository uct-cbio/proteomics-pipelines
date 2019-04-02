#!/usr/bin/env bash

set -e

sudo ./build.sh

outfolder=/data/bio/mtb/proteogenomics/17_3_missed
interproscan_path="/data/bio/interproscan/5_31_70/interproscan-5.31-70.0"

sudo docker run -it --rm -v ${outfolder}:/root/output \
    -v $HOME/repos/proteomics-pipelines/bin/config/mq_proteogeomics_example_config_2strains.yml:/root/config.yaml \
    -v /data/bio/mtb/gbrowse/data/gff3/Mycobacterium_tuberculosis_h37rv.ASM19595v2.34.chromosome.Chromosome.gff3:/root/features.gff3 \
    -v /data/bio/mtb/mq/17_3_missed/txt:/root/txt \
    -v /data/bio/mtb/genomes:/root/genomes \
    -v /data/bio/mtb/operons/3354.opr:/root/3354.txt \
    -v ${interproscan_path}:/interproscan \
    cbio/mqproteogenomics:latest

# Load data to elasticsearch
#cd ../csv2es && sudo ./run.sh localhost 9230 "s507_s5527" ${outfolder}

