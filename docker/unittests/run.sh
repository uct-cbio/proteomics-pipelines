#!/usr/bin/env bash
sudo ./build.sh
sudo docker run -it --rm -v $HOME/bio:/root/bio -v $HOME/repos/proteomics-pipelines:/root/proteomics-pipelines -v $HOME/repos/cbio-proteogenomics-tests:/root/unittests cbio/unittests:latest 
