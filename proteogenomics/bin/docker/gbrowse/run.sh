#!/usr/bin/env bash

./build.sh
docker run -it --rm --net='host' cbio/gbrowse2:latest
