#/usr/bin/env bash
set -e

# Change build context to the root of the directory
version=v3.0.6

name=bio
docker pull thyscbio/${name}:${version} && echo "version exists" && exit 1 || echo "Gonna build."
cd ../..
docker build --build-arg ftp_proxy=$ftp_proxy --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy -f docker/${name}/Dockerfile  -t thyscbio/${name}:${version} .
docker login https://index.docker.io/v1/thyscbio
docker push  thyscbio/${name}:${version} || exit 1

docker tag thyscbio/${name}:${version} thyscbio/${name}:latest
docker push  thyscbio/${name}:latest
