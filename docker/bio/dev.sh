#!/bin/bash
set -a
version=v2.3
home=/home/$(who -m | awk '{print $1;}')
config=${home}/.proteomics-pipelines_profile
source $config
docker run -it --rm -v ${CBIO_PROTEOGENOMICS_TESTS_PATH}:${CBIO_PROTEOGENOMICS_TESTS_PATH} -v ${PROTEOMICS_PIPELINES_PATH}:${PROTEOMICS_PIPELINES_PATH} thyscbio/bio:${version}  /bin/bash 
