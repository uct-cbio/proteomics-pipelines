#!/usr/bin/env bash

./build.sh
docker run -it --rm -v /data/bio/interproscan/5_31_70/interproscan-5.31-70.0:/interproscan \
    --net='host' cbio/iprscan:latest /bin/bash
