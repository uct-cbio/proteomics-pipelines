#!/usr/bin/env bash
sudo ./build.sh
sudo docker run -it --rm -v $HOME/bio/mtb/proteogenomics/multistrain1:/root/output \
    -v $HOME/repos/cbio-pipelines/proteogenomics/bin/config/mq_proteogeomics_example_config.yml:/root/config.yaml \
    -v $HOME/bio/gbrowse/data/gff3/Mycobacterium_tuberculosis_h37rv.ASM19595v2.34.chromosome.Chromosome.gff3:/root/features.gff3 \
    -v $HOME/bio/mtb/mq/multistrain1/txt:/root/txt \
    -v $HOME/bio/mtb/genomes:/root/genomes \
    cbio/mqproteogenomics:latest 
