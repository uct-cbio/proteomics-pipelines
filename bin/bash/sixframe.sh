#!/usr/bin/env bash

set -e

input=$1
config=$2

mkdir -p ${input}/sixframe 



cp  ${config} ${input}/sixframe/

for filename in ${input}/*.fasta; do
    base=${input}/sixframe/$(basename ${filename%.fasta})
    sf.py ${filename} $(basename ${filename%.fasta}) ${config} ${base}
    cat ${input}/sixframe/*_sixframe_prot.fasta > ${input}/sixframe/combined_prot.fasta
    cat ${input}/sixframe/*_sixframe_nuc.fasta > ${input}/sixframe/combined_nuc.fasta
done

