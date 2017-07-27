#!/usr/bin/env bash

txt=$1
output=${txt}/cbio_metaproteomics
echo ${txt}
#rm -rf $output

mkdir $output

mqmetaproteomics.py ${txt} ${output}
