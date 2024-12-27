#!/bin/bash
set -a
version=v3.0.19

home=/home/$(who -m | awk '{print $1;}')
config=${home}/.proteomics-pipelines_profile
source $config

docker run  -it --rm -p 8181:9090 -v ${CBIO_PROTEOGENOMICS_TESTS_PATH}:${CBIO_PROTEOGENOMICS_TESTS_PATH} -v ${PROTEOMICS_PIPELINES_PATH}/bin:/home/bio/bin -v ${PROTEOMICS_PIPELINES_PATH}/lib:/home/bio/lib thyscbio/bio:${version}  /bin/bash 
