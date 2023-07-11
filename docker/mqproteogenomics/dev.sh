#!/bin/bash
set -a
version=v1.22.2
home=/home/$(who -m | awk '{print $1;}')
config=${home}/.proteomics-pipelines_profile
source $config
docker run -it --rm -p 8181:9090 -v ${CBIO_PROTEOGENOMICS_TESTS_PATH}:/home/testing/cbio-proteogenomics-tests -v ${PROTEOMICS_PIPELINES_PATH}:/home/testing/proteomics-pipelines thyscbio/mqproteogenomics:${version}  /bin/bash 
