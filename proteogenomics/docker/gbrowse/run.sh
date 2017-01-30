#!/usr/bin/env bash

./build.sh
#docker run -it --rm --net='host' -v $HOME/bio/gbrowse/data:/data cbio/gbrowse2:latest
docker run -d --restart "always" -p 8080:80 -v $HOME/bio/gbrowse/data:/data --name cbio_gbrowse2 cbio/gbrowse2:latest 
