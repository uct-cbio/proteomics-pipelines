#!/usr/bin/env bash

input=$1
output=$2

set -e

mkdir $2
cp $1/* $2

cd $2 && tar -zxvf * && cat *.fasta > uniprot_current.fasta 
