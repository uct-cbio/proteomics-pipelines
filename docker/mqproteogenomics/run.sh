#!/usr/bin/env bash
set -e

sudo ./build.sh

outfolder=$HOME/bio/mtb/proteogenomics/17_3_missed

sudo docker run -it --rm -v ${outfolder}:/root/output \
    -v $HOME/repos/cbio-pipelines/proteomics/bin/config/mq_proteogeomics_example_config_2strains.yml:/root/config.yaml \
    -v $HOME/bio/gbrowse/data/gff3/Mycobacterium_tuberculosis_h37rv.ASM19595v2.34.chromosome.Chromosome.gff3:/root/features.gff3 \
    -v $HOME/bio/mtb/mq/17_3_missed/txt:/root/txt \
    -v $HOME/bio/mtb/all_genomes:/root/genomes \
    -v $HOME/bio/mtb/operons/3354.opr:/root/3354.txt \
    cbio/mqproteogenomics:latest

# Load data to elasticsearch
#cd ../csv2es && sudo ./run.sh localhost 9230 "s507_s5527" ${outfolder}

