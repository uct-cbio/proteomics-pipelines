#!/usr/bin/env bash

ES_HOST=$1
ES_PORT=$2
ES_ALIAS=$3
ES_USER=$4
ES_PASSWORD=$5
ES_DATA=$6
./build.sh
sudo docker run -it --rm -v $ES_DATA:/root/data -e ES_HOST=$ES_HOST \
    -e ES_PORT=$ES_PORT \
    -e ES_ALIAS=$ES_ALIAS \
    -e ES_USER=$ES_USER \
    -e ES_PASSWORD=$ES_PASSWORD \
    cbio/csv2es:latest
