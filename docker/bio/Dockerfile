FROM amd64/r-base:latest

ENV http_proxy ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

RUN apt-get update && apt install -y default-jdk

#RUN wget https://launchpad.net/~openjdk-r/+archive/ubuntu/ppa/+build/8156693/+files/openjdk-8-jdk_8u66-b17-1_amd64.deb
#RUN dpkg -i openjdk-8-jre_8u66-b17-1_amd64.deb || apt-get -y install -f && dpkg -i openjdk-8-jdk_8u66-b17-1_amd64.deb
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install software-properties-common

RUN java -version
RUN apt-get update && apt-get -y  install dirmngr

# PYTHON 3
RUN apt-get install -y python3-dev python3-pip
RUN pip3 install --upgrade pip
RUN pip3 install pandas scipy numpy virtualenv
RUN pip3 install BioPython rpy2

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get install -y curl libssl-dev libxml2-dev libcurl4-openssl-dev libssh2-1-dev openssl curl openssl
RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('Rcpp', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('mime', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('optparse', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('cluster',dependencies=TRUE)"
RUN apt-get update && apt-get install -f -y libglvnd-dev libglx-dev  libgl-dev  libgl1-mesa-dev octave liboctave-dev
RUN apt-get install -y openmpi-bin openmpi-doc libopenmpi-dev
RUN R -e "install.packages('Rmpi', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "BiocManager::install( c('NMF', 'Matrix', 'irlba'), dependencies=TRUE)"

RUN apt-get install -y mesa-common-dev libglu1-mesa-dev freeglut3-dev bwidget libfreetype6-dev
RUN R -e "BiocManager::install( c('igraph', 'lintr', 'curl'), dependencies=TRUE)"

RUN apt-get install -y libnetcdf-dev  
RUN R -e "BiocManager::install( c('mzR', 'kernlab', 'robustbase', 'mclust'), dependencies=TRUE)"

RUN apt-get -y install ed
RUN R -e "BiocManager::install( c('MLInterfaces', 'S4Vectors', 'htmltools', 'preprocessCore', 'pRoloc', 'MSmbase', 'rTANDEM'), dependencies=TRUE)"

RUN R CMD javareconf
RUN R -e "BiocManager::install( c('rJava', 'RforProteomics', 'gage', 'limma', 'ggplot', 'dendextend', 'mvtnorm', 'tmvtnorm'), dependencies=TRUE)"
RUN R -e "library(devtools); install_github('vqv/ggbiplot')"

RUN R -e "BiocManager::install( c('imputeLCMD', 'pvclust', 'colorspace', 'data.table', 'qvalue', 'dplyr'), dependencies=TRUE)"

RUN apt-get install -y libcairo2-dev
RUN R -e "BiocManager::install( c('gdtools', 'rngtools', 'ggplot2'), dependencies=TRUE)"

RUN pip3 install matplotlib
RUN R -e "install.packages('epitools', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('plotrix', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('sciplot', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('dunn.test', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('FSA',,'http://www.rforge.net/', dependencies=TRUE) "
RUN pip3 install urllib3

# Ruby and unipept
RUN apt-get update && apt-get -y install ruby-full
RUN gem install unipept

RUN groupadd -r bio && useradd -r -g bio bio
WORKDIR /home/bio/

# Install BLAST
RUN apt-get -y install wget
ENV blast_ftp ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz
RUN mkdir /home/bio/blast && cd /home/bio/blast && wget ${blast_ftp} && tar -zxvf * && rm -rf *.tar.gz 
RUN cd /home/bio/blast/ncbi*/bin && cp * /usr/local/bin

# Clustalw
ENV clustalw_ftp ftp://ftp.ebi.ac.uk/pub/software/clustalw2/2.1/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
RUN mkdir /home/bio/clustalw && cd /home/bio/clustalw && wget ${clustalw_ftp} && tar -zxvf * && rm -rf *.tar.gz
RUN cd /home/bio/clustalw/clustalw* && cp * /usr/local/bin

# MUSCLE
ENV muscle_ftp http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux32.tar.gz
RUN mkdir /home/bio/muscle && cd /home/bio/muscle && wget ${muscle_ftp} && tar -zxvf * && rm -rf *.tar.gz
RUN cd /home/bio/muscle/ && cp muscle* /usr/local/bin/muscle

# Python deps
RUN pip3 install natsort 
RUN pip3 install PyYAML
RUN apt-get update
RUN pip3 install matplotlib_venn
RUN apt-get update
RUN apt-get -y install python3-tk python3-cffi

RUN apt-get -y install libffi-dev
RUN pip3 install --upgrade cffi
RUN pip3 install --upgrade pip
RUN pip3 install rpy2
RUN pip3 install datrie

RUN pip3 install scikit-posthocs
RUN R -e "install.packages('Rcpp')"
#RUN R -e "devtools::install_version('car', version = '2.1-6', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('car', dependencies=TRUE, repos='http://cran.us.r-project.org')"

#RUN R -e "install.packages('tvtnorm', type='source', repos=NULL) "
RUN R -e "install.packages('FSA') "

RUN apt-get -y install libdw-dev

RUN pip3 install mygene xmltramp2
# Issue with openblas
RUN apt-get -y install git parallel

RUN git clone https://github.com/bmbolstad/preprocessCore.git && cd preprocessCore/ && R CMD INSTALL --configure-args="--disable-threading"  .

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# Leave to end
#ENV PATH $PATH::/home/bio/bin/python:/home/bio/bin/R:/home/bio/bin/bash:/home/bio/bin/python
#ENV PYTHONPATH $PYTHONPATH::/home/bio/lib

#COPY lib /home/bio/lib
#COPY bin /home/bio/bin


