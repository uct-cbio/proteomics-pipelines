#!/usr/bin/env bash
sudo ./build.sh
sudo docker run -it --rm -v /home/thys/bio/mtb/genomes:/root/fasta  cbio/sixframe:latest 
