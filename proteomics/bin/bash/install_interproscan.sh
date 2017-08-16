#!/usr/bin/env bash

set -e
mkdir $1
cd $1

wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.19-58.0/interproscan-5.19-58.0-64-bit.tar.gz && wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.19-58.0/interproscan-5.19-58.0-64-bit.tar.gz.md5 && md5sum -c interproscan-5.19-58.0-64-bit.tar.gz.md5
# Must return *interproscan-5.19-58.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.
tar -pxvzf interproscan-5.19-58.0-*-bit.tar.gz

cd interproscan-5.19-58.0/data && wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-10.0.tar.gz && wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-10.0.tar.gz.md5 && md5sum -c panther-data-10.0.tar.gz.md5
# This must return *panther-data-10.0.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.
tar -pxvzf panther-data-10.0.tar.gz






