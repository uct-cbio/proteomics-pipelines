#!/bin/bash
set -a
version=v3.0.5
home=/home/$(who -m | awk '{print $1;}')
config=${home}/.proteomics-pipelines_profile
source $config

DEV_PYTHONPATH=${PROTEOMICS_PIPELINES_PATH}/lib/

DEV_PATH="$PATH::${PROTEOMICS_PIPELINES_PATH}/bin/python/:${PROTEOMICS_PIPELINES_PATH}/bin/R:${PROTEOMICS_PIPELINES_PATH}/bin/bash/:${PROTEOMICS_PIPELINES_PATH}/bin/perl/"

docker run -it --rm -p 8181:9090 -e PYTHONPATH=${DEV_PYTHONPATH} -e PATH=${DEV_PATH} -v ${CBIO_PROTEOGENOMICS_TESTS_PATH}:${CBIO_PROTEOGENOMICS_TESTS_PATH} -v ${PROTEOMICS_PIPELINES_PATH}:${PROTEOMICS_PIPELINES_PATH} thyscbio/bio:${version}  /bin/bash 
