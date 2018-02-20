#!/usr/bin/env bash
set -e

sudo ./build.sh

outfolder=$HOME/bio/mtb/proteogenomics/multistrain1_test

sudo docker run -it --rm -v ${outfolder}:/root/output \
    -v $HOME/repos/cbio-pipelines/proteogenomics/bin/config/mq_proteogeomics_example_config.yml:/root/config.yaml \
    -v $HOME/bio/mtb/mq/multistrain1_test/txt:/root/txt \
    -v $HOME/bio/gbrowse/data/gff3/Mycobacterium_tuberculosis_h37rv.ASM19595v2.34.chromosome.Chromosome.gff3:/root/features.gff3 \
    -v $HOME/bio/mtb/genomes:/root/genomes \
    -v $HOME/bio/mtb/operons/3354.opr:/root/3354.txt \
    cbio/mqproteogenomics:latest

# Load data to elasticsearch
cd ../csv2es && sudo ./run.sh 192.168.7.4 9210 "test" "elastic" "changeme" ${outfolder}
 

