FROM ubuntu:latest


ENV  DEBIAN_FRONTEND noninteractive 
RUN apt-get update && apt-get -y install software-properties-common default-jdk dirmngr wget libdw-dev git parallel 
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh && bash Miniconda*.sh -b && rm Miniconda*.sh
ENV PATH=/root/miniconda3/bin:${PATH} 
RUN conda update -y conda
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --add channels r
RUN conda install pandas scipy numpy 
RUN conda install virtualenv xmltramp2 cffi 
RUN conda install rpy2 datrie natsort 
RUN conda install PyYAML urllib3 scikit-posthocs matplotlib
RUN conda install -c bioconda blast muscle clustalw BioPython mygene
RUN conda install -c conda-forge matplotlib-venn
#RUN conda install -c r r r-essentials


RUN apt-get update && apt-get install -y curl libssl-dev libxml2-dev libcurl4-openssl-dev libssh2-1-dev openssl curl openssl
RUN apt-get update && apt-get install -f -y libglvnd-dev libglx-dev  libgl-dev  libgl1-mesa-dev octave liboctave-dev
RUN apt-get update && apt-get install -y mesa-common-dev libglu1-mesa-dev freeglut3-dev bwidget libfreetype6-dev
RUN apt-get update && apt-get install -y libnetcdf-dev  libboost-all-dev
RUN apt-get update && apt-get install -y gcc gfortran libnetcdf-dev libnetcdff-dev libhdf5-serial-dev
RUN apt-get update && apt-get install -y openmpi-bin openmpi-doc libopenmpi-dev
RUN apt-get update && apt-get -y install zlib1g
RUN apt-get update && apt-get install -y libcairo2-dev
RUN apt-get update && apt-get -y install ed
RUN apt-get update && apt-get -y install libffi-dev libxml2-dev
RUN conda install -vv -c bioconda bioconductor-mzr

RUN conda install -c conda-forge r-rcpp r-mime r-optparse r-devtools r-rmpi
RUN conda install -c r r-cluster
RUN conda install -c conda-forge r-nmf r-matrix r-irlba r-igraph r-lintr r-curl
RUN conda install -c bioconda bioconductor-rhdf5lib
RUN conda install -c conda-forge r-kernlab r-robustbase r-mclust
#RUN R -e "BiocManager::install('Rhdf5lib')"
#RUN R -e "BiocManager::install('mzR')"
#RUN R -e "BiocManager::install( c('mzR', 'kernlab', 'robustbase', 'mclust'), dependencies=TRUE)"
RUN conda install -c bioconda bioconductor-mlinterfaces
RUN conda install -c bioconda bioconductor-s4vectors
RUN conda install -c conda-forge r-htmltools

RUN conda install -c bioconda bioconductor-preprocesscore
RUN conda install -c bioconda bioconductor-proloc
RUN conda install -c bioconda bioconductor-msnbase
#RUN conda install -c bioconda bioconductor-rtandem

#RUN R -e "BiocManager::install( c('MLInterfaces', 'S4Vectors', 'htmltools', 'preprocessCore', 'pRoloc', 'MSmbase', 'rTANDEM'), dependencies=TRUE)"

RUN R CMD javareconf
RUN conda install -c r r-rjava
#RUN R -e "BiocManager::install( 'RforProteomics', dependencies=TRUE)"
#RUN conda install -c bioconda bioconductor-rforproteomics
RUN conda install -c bioconda bioconductor-gage
RUN conda install -c bioconda bioconductor-limma
RUN conda install -c conda-forge r-ggplot2
RUN conda install -c conda-forge r-dendextend r-mvtnorm r-tmvtnorm

RUN conda install -c bioconda r-ggbiplot
#RUN R -e "library(devtools); install_github('vqv/ggbiplot')"
RUN conda install -c bioconda r-imputelcmd
RUN conda install -c conda-forge r-pvclust
RUN conda install -c conda-forge r-colorspace
RUN conda install -c conda-forge r-data.table
RUN conda install -c bioconda bioconductor-qvalue
RUN conda install -c conda-forge r-dplyr
RUN conda install -c conda-forge r-gdtools
RUN conda install -c conda-forge r-rngtools
#RUN R -e "BiocManager::install( c('imputeLCMD', 'pvclust', 'colorspace', 'data.table', 'qvalue', 'dplyr'), dependencies=TRUE)"

#RUN R -e "BiocManager::install( c('gdtools', 'rngtools', 'ggplot2'), dependencies=TRUE)"
#RUN conda install -c bird r-epitools
RUN conda install -c conda-forge r-plotrix
RUN conda install -c conda-forge r-sciplot
RUN conda install -c bioconda r-dunn.test
RUN conda install -c conda-forge r-fsa
RUN conda install -c bioconda r-car
RUN conda install -c r r-gplots
RUN conda install -c bioconda r-corrplot
RUN conda install -c bioconda bioconductor-pathview
##RUN R -e "install.packages('epitools', dependencies=TRUE, repos='http://cran.us.r-project.org')"
#RUN R -e "install.packages('plotrix', dependencies=TRUE, repos='http://cran.us.r-project.org')"
#RUN R -e "install.packages('sciplot', dependencies=TRUE, repos='http://cran.us.r-project.org')"
#RUN R -e "install.packages('dunn.test', dependencies=TRUE, repos='http://cran.us.r-project.org')"
#RUN R -e "install.packages('FSA',,'http://www.rforge.net/', dependencies=TRUE) "

# Ruby and unipept
RUN apt-get update && apt-get -y install ruby-full
RUN gem install unipept

RUN groupadd -r bio && useradd -r -g bio bio
WORKDIR /home/bio/

# Install BLAST
#ENV blast_ftp ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz
#RUN mkdir /home/bio/blast && cd /home/bio/blast && wget ${blast_ftp} && tar -zxvf * && rm -rf *.tar.gz 
#RUN cd /home/bio/blast/ncbi*/bin && cp * /usr/local/bin

# Clustalw
#ENV clustalw_ftp ftp://ftp.ebi.ac.uk/pub/software/clustalw2/2.1/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
#RUN mkdir /home/bio/clustalw && cd /home/bio/clustalw && wget ${clustalw_ftp} && tar -zxvf * && rm -rf *.tar.gz
#RUN cd /home/bio/clustalw/clustalw* && cp * /usr/local/bin

# MUSCLE
#ENV muscle_ftp http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux32.tar.gz
#RUN mkdir /home/bio/muscle && cd /home/bio/muscle && wget ${muscle_ftp} && tar -zxvf * && rm -rf *.tar.gz
#RUN cd /home/bio/muscle/ && cp muscle* /usr/local/bin/muscle

# Python deps
#RUN apt-get -y install python3-tk python3-cffi

#RUN R -e "install.packages('Rcpp')"
#RUN R -e "devtools::install_version('car', version = '2.1-6', repos = 'http://cran.us.r-project.org')"
#RUN R -e "devtools::install_version('car', dependencies=TRUE, repos='http://cran.us.r-project.org')"

#RUN R -e "install.packages('tvtnorm', type='source', repos=NULL) "
#RUN R -e "install.packages('FSA') "


RUN git clone https://github.com/bmbolstad/preprocessCore.git && cd preprocessCore/ && R CMD INSTALL --configure-args="--disable-threading"  .
#RUN rm /bin/sh && ln -s /bin/bash /bin/sh

RUN apt-get update && apt-get install -y libssl-dev openssl
RUN conda update conda && conda install -c bioconda bcftools minimap2  
RUN conda config --add channels bioconda && conda config --add channels conda-forge 
RUN conda install samtools  
RUN conda install tabix  
RUN conda install jbrowse2  
RUN conda install emboss  
RUN yes | conda update --all  
RUN conda init 
RUN apt-get update && apt-get install -y software-properties-common build-essential zlib1g-dev  genometools 


#ENV PYTHONPATH /root/miniconda3/lib/python3.9/site-packages::$PYTHONPATH
# testing paths
#ENV PYTHONPATH $PYTHONPATH::/home/testing/proteomics-pipelines/lib
# production paths
ENV PYTHONPATH $PYTHONPATH::/home/bio/lib

# testing paths (if files are mounted it will use those
#ENV PATH $PATH::/home/testing/proteomics-pipelines/bin/python:/home/testing/proteomics-pipelines/bin/R:/home/testing/proteomics-pipelines/bin/bash:/home/testing/proteomics-pipelines/bin/python
# production paths
ENV PATH $PATH::/home/bio/bin/python:/home/bio/bin/R:/home/bio/bin/bash:/home/bio/bin/python


#COPY lib /home/bio/lib
#COPY bin /home/bio/bin


