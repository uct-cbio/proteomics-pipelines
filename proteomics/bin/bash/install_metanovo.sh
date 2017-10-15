#!/usr/bin/env bash

peptideshaker='http://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/1.16.2/PeptideShaker-1.16.2.zip'

searchgui='http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.2.20/SearchGUI-3.2.20-mac_and_linux.tar.gz'

denovogui='http://genesis.ugent.be/maven2/com/compomics/denovogui/DeNovoGUI/1.15.11/DeNovoGUI-1.15.11-mac_and_linux.tar.gz'

utilities='http://genesis.ugent.be/maven2/com/compomics/utilities/4.11.19/utilities-4.11.19.zip'

proteoannotator='http://www.proteoannotator.org/datasets/releases/ProteoAnnotator-1.7.86.zip'

ruby='https://cache.ruby-lang.org/pub/ruby/2.2/ruby-2.2.5.tar.gz'

proteowizard='https://sourceforge.net/projects/proteowizard/files/latest/download'

blast='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.5.0+-x64-linux.tar.gz'

clustalw='http://www.clustal.org/download/current/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz'

muscle='http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz'

kibana='https://download.elastic.co/kibana/kibana/kibana-4.5.4-linux-x64.tar.gz'

elasticsearch='https://download.elastic.co/elasticsearch/release/org/elasticsearch/distribution/tar/elasticsearch/2.3.5/elasticsearch-2.3.5.tar.gz'

R='https://cran.r-project.org/src/base/R-3/R-3.2.3.tar.gz'

cmake='https://cmake.org/files/v3.7/cmake-3.7.0-rc2.tar.gz'

llvm='http://llvm.org/releases/3.7.0/llvm-3.7.0.src.tar.xz'
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

#############
## CLUSTALW # 
#############

url=$clustalw
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

#############
## MUSCLE   #
#############
url=$muscle
file="${url##*/}"
name="${file%.tar.gz}"
base=muscle
if [ ! -d "$dir/$base" ]; then 
    mkdir $dir/$base
fi 
cd $dir/$base 
echo $name
if [ ! -d "$name" ]; then 
    mkdir $name && cd $name && wget $url && tar -zxvf *.tar.gz && rm -rf *.tar.gz && cp * muscle 
fi 

##################
## Elasticsearch #
##################
url=$elasticsearch
file="${url##*/}"
name="${file%.tar.gz}"
base=elasticsearch
if [ ! -d "$dir/$base" ]; then 
    mkdir $dir/$base
fi 
cd $dir/$base 
echo $name
if [ ! -d "$name" ]; then 
    mkdir $name && cd $name && wget $url && tar -zxvf *.tar.gz && rm -rf *.tar.gz  
fi 


#################
## Kibana       #
#################
url=$kibana
file="${url##*/}"
name="${file%.tar.gz}"
base=kibana
if [ ! -d "$dir/$base" ]; then 
    mkdir $dir/$base
fi 
cd $dir/$base 
echo $name
if [ ! -d "$name" ]; then 
    mkdir $name && cd $name && wget $url && tar -zxvf *.tar.gz && rm -rf *.tar.gz  
fi 


#################
## R ###########
#################
url=$R
file="${url##*/}"
name="${file%.tar.gz}"
base=R
if [ ! -d "$dir/$base" ]; then 
    mkdir $dir/$base
fi 
cd $dir/$base 
echo $name
if [ ! -d "$name" ]; then 
    mkdir $name && cd $name && wget $url && tar -zxvf *.tar.gz && rm -rf *.tar.gz && cd R* && ./configure --prefix=$dir/$base/$name && make && make install 
fi 

##############
# cmake      #
##############
url=$cmake
file="${url##*/}"
name="${file%.tar.gz}"
base=cmake
if [ ! -d "$dir/$base" ]; then 
    mkdir $dir/$base
fi 
cd $dir/$base 
echo $name
if [ ! -d "$name" ]; then 
    mkdir $name && cd $name && wget $url && tar -xvf *.tar.* && rm -rf *.tar.* && cd $name* && ./configure --prefix=$dir/$base/$name && make && make install 
fi 

##############
# llvm       #
##############
url=$llvm
file="${url##*/}"
name="${file%.src.tar.xz}"
base=llvm
if [ ! -d "$dir/$base" ]; then 
    mkdir $dir/$base
fi 
cd $dir/$base 
echo $name
if [ ! -d "$name" ]; then 
    rm -rf build && mkdir build
    rm -rf dest && mkdir dest
    mkdir $name && cd $name && wget $url && tar -xvf *.tar.xz && rm -rf *.tar.xz && cd $dir/$base/build && cmake $dir/$base/$name/$name* && cmake -DCMAKE_INSTALL_PREFIX:PATH=$dir/$base/dest . && make all install # cmake $dir/$base/$name/$name* && cmake --build . --target install #|| rm -rf $dir/$base/$name
fi 














