#!/usr/bin/env bash

./build.sh
docker run -it --rm --net='host' --name cbio_gbrowse2 -v $HOME/bio/gbrowse/data:/data cbio/gbrowse2:latest
#docker run -it --rm --net='host' -v $HOME/bio/gbrowse/data:/data cbio/gbrowse2:latest /bin/bash
