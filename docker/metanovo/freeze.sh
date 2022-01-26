#/usr/bin/env bash


metanovo_version=v1.8

sudo docker run -it --rm thyscbio/metanovo:${metanovo_version} pip freeze
