#/usr/bin/env bash

# Change build context to the root of the directory

cd ../..

sudo docker login https://index.docker.io/v1/thyscbio

sudo docker build --build-arg ftp_proxy=$ftp_proxy --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy -f docker/mqmetaproteomics/Dockerfile  -t thyscbio/mqmetaproteomics:latest .
sudo docker push  thyscbio/mqmetaproteomics:latest
