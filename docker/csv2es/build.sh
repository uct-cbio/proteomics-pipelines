#!/usr/bin/env bash

# Change build context
cd ../..

sudo docker build --build-arg ftp_proxy=$ftp_proxy --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy -f docker/csv2es/Dockerfile -t cbio/csv2es:latest . 
