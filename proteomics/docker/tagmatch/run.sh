#/usr/bin/env bash

MGF_FOLDER=''
CONFIG_FILE=''

./build.sh
sudo docker run cbio/tagmatch:latest
