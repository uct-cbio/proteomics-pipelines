#!/usr/bin/env bash

searchgui='http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.2.20/SearchGUI-3.2.20-mac_and_linux.tar.gz'

denovogui='http://genesis.ugent.be/maven2/com/compomics/denovogui/DeNovoGUI/1.15.11/DeNovoGUI-1.15.11-mac_and_linux.tar.gz'

utilities='http://genesis.ugent.be/maven2/com/compomics/utilities/4.12.0/utilities-4.12.0.zip'


#blast='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.5.0+-x64-linux.tar.gz'

#############################
# Create the base directory #
#############################

set -e
dir=${METANOVO_DEPENDENCIES}
echo "Installing dependencies to ${dir}"
if [ ! -d "$dir" ]; then 
    mkdir $dir
fi 


#################
# utilties      #
#################

url=$utilities
file="${url##*/}"
name="${file%.zip}"
base=$(echo $name | cut -f 1 -d '-')
if [ ! -d "$dir/$base" ]; then 
    mkdir $dir/$base
fi 
cd $dir/$base 
echo $name
if [ ! -d "$name" ]; then 
    wget $url && unzip *.zip && rm -rf *.zip 
fi 


#############
# SearchGUI #
#############

url=$searchgui
file="${url##*/}"
name="${file%-mac_and_linux.tar.gz}"
base=$(echo $name | cut -f 1 -d '-')
if [ ! -d "$dir/$base" ]; then 
    mkdir $dir/$base
fi 
cd $dir/$base 
echo $name
if [ ! -d "$name" ]; then 
    wget $url && tar -zxvf *.tar.gz && rm -rf *.tar.gz
fi 
#############
# DeNovoGUI #
#############

url=$denovogui
file="${url##*/}"
name="${file%-mac_and_linux.tar.gz}"
base=$(echo $name | cut -f 1 -d '-')
if [ ! -d "$dir/$base" ]; then 
    mkdir $dir/$base
fi 
cd $dir/$base 
echo $name
if [ ! -d "$name" ]; then 
    wget $url && tar -zxvf *.tar.gz && rm -rf *.tar.gz
fi 






