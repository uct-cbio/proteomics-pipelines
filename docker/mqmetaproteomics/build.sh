#/usr/bin/env bash
set -e

# Change build context to the root of the directory
version=v2.0

name=mqmetaproteomics
docker pull thyscbio/${name}:${version} && echo "version exists" && exit 1 || echo "Gonna build."
cd ../..
sudo docker build --build-arg ftp_proxy=$ftp_proxy --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy -f docker/mqmetaproteomics/Dockerfile  -t thyscbio/${name}:${version} .
sudo docker login https://index.docker.io/v1/thyscbio
sudo docker push  thyscbio/${name}:${version}
