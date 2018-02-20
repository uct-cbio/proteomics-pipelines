#!/usr/bin/env bash

peptideshaker='http://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/1.16.2/PeptideShaker-1.16.2.zip'

searchgui='http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.2.20/SearchGUI-3.2.20-mac_and_linux.tar.gz'

denovogui='http://genesis.ugent.be/maven2/com/compomics/denovogui/DeNovoGUI/1.15.11/DeNovoGUI-1.15.11-mac_and_linux.tar.gz'

utilities='http://genesis.ugent.be/maven2/com/compomics/utilities/4.11.19/utilities-4.11.19.zip'

proteoannotator='http://www.proteoannotator.org/datasets/releases/ProteoAnnotator-1.7.86.zip'

ruby='https://cache.ruby-lang.org/pub/ruby/2.2/ruby-2.2.5.tar.gz'

proteowizard='https://sourceforge.net/projects/proteowizard/files/latest/download'

blast='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.5.0+-x64-linux.tar.gz'

#############################
# Create the base directory #
#############################

set -e
dir=$HOME/software
if [ ! -d "$dir" ]; then 
    mkdir $dir
fi 

#################
# PeptideShaker #
#################

url=$peptideshaker
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

###################
# proteoannotator #
###################

url=$proteoannotator
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

########
# ruby #
########

url=$ruby
file="${url##*/}"
name="${file%.tar.gz}"
base=$(echo $name | cut -f 1 -d '-')
if [ ! -d "$dir/$base" ]; then 
    mkdir $dir/$base
fi 
cd $dir/$base 
echo $name
if [ ! -d "$name" ]; then 
    mkdir $name && cd $name && wget $url && tar -zxvf *.tar.gz && rm -rf *.tar.gz && cd ruby* && ./configure --prefix=$dir/$base/$name && make && make install 
fi 

##########
## BLAST #
##########

url=$blast
file="${url##*/}"
name="${file%.tar.gz}"
base=$(echo $name | cut -f 1 -d '-')
if [ ! -d "$dir/$base" ]; then 
    mkdir $dir/$base
fi 
cd $dir/$base 
echo $name
if [ ! -d "$name" ]; then 
    mkdir $name && cd $name && wget $url && tar -zxvf *.tar.gz && rm -rf *.tar.gz  
fi 





