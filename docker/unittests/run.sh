#!/usr/bin/env bash
sudo ./build.sh
sudo docker run -it --rm -v $HOME/bio:/root/bio -v $HOME/repos/cbio_proteomics:/root/cbio_proteomics cbio/unittests:latest 
