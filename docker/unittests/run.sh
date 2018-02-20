#!/usr/bin/env bash
sudo ./build.sh
sudo docker run -it --rm -v $HOME/bio:/root/bio -v $HOME/repos/cbio-proteogenomics-tests:/root/unittests -v $HOME/repos/cbio-pipelines:/root/cbio-pipelines cbio/unittests:latest 
