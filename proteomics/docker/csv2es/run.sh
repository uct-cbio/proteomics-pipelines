#!/usr/bin/env bash

ES_HOST=$1
ES_PORT=$2
ES_ALIAS=$3
ES_DATA=$4
./build.sh
sudo docker run -it --rm --net=host -v $ES_DATA:/root/data -e ES_HOST=$ES_HOST \
    -e ES_PORT=$ES_PORT \
    -e ES_ALIAS=$ES_ALIAS \
    cbio/csv2es:latest
