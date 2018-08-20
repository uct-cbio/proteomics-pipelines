FROM ubuntu:xenial

ENV http_proxy ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

RUN apt-get update

# PYTHON 3
RUN apt-get install -y python3-dev python3-pip
RUN pip3 install --upgrade pip
RUN pip3 install pandas scipy numpy virtualenv
RUN pip3 install BioPython

# PYTHON 2.7
RUN apt-get install -y python-dev python-pip
RUN virtualenv -p python2.7 /root/ve27

# Tools
RUN apt-get install -y software-properties-common

# Install Oracle Java 8
RUN \
  echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
  add-apt-repository -y ppa:webupd8team/java && \
  apt-get update && \
  apt-get install -y oracle-java8-installer && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /var/cache/oracle-jdk8-installer
ENV JAVA_HOME /usr/lib/jvm/java-8-oracle

# R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb [arch=amd64,i386] http://cran.rstudio.com/bin/linux/ubuntu xenial/'
RUN apt-get update && apt-get install -y r-base r-base-dev
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite()"
RUN R -e "install.packages('Rcpp', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('mime', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('optparse', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN apt-get install -y curl libssl-dev libxml2-dev libcurl4-openssl-dev libssh2-1-dev 
RUN R -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('cluster',dependencies=TRUE)"
RUN apt-get install -y octave liboctave-dev
RUN apt-get install -y r-cran-rmpi openmpi-bin openmpi-doc libopenmpi-dev
RUN R -e "install.packages('Rmpi', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('NMF',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('Matrix',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('irlba',dependencies=TRUE)"
RUN apt-get install -y mesa-common-dev libglu1-mesa-dev freeglut3-dev bwidget
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('igraph',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('lintr',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('curl',dependencies=TRUE)"
RUN apt-get install -y libnetcdf-dev libnetcdf11  
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('mzR',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('kernlab',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('robustbase',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('mclust',dependencies=TRUE)"
RUN apt-get -y install ed
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('MLInterfaces',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('S4Vectors',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('htmltools',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('preprocessCore',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('pRoloc',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('MSnbase',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('rTANDEM',dependencies=TRUE)"
RUN R CMD javareconf
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('rJava',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('RforProteomics',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('gage',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('limma',dependencies=TRUE)"
RUN R -e "library(devtools); install_github('vqv/ggbiplot')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('gplots',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('dendextend',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('imputeLCMD',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('pvclust',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('colorspace',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('data.table',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('qvalue',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('dplyr',dependencies=TRUE)"
RUN apt-get install -y libcairo2-dev
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('gdtools',dependencies=TRUE)"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('ggplot2',dependencies=TRUE)"
RUN pip3 install matplotlib
RUN pip3 install rpy2
RUN R -e "install.packages('epitools', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('plotrix', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('sciplot', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('dunn.test', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('FSA',,'http://www.rforge.net/', dependencies=TRUE) "
RUN pip3 install urllib3

# Ruby and unipept
RUN apt-get update && apt-get -y install ruby-full
RUN gem install unipept

# Install BLAST
ENV blast_ftp ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-x64-linux.tar.gz
RUN mkdir /root/blast && cd /root/blast && wget ${blast_ftp} && tar -zxvf * && rm -rf *.tar.gz 
RUN cd /root/blast/ncbi*/bin && cp * /usr/local/bin

# Clustalw
ENV clustalw_ftp ftp://ftp.ebi.ac.uk/pub/software/clustalw2/2.1/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
RUN mkdir /root/clustalw && cd /root/clustalw && wget ${clustalw_ftp} && tar -zxvf * && rm -rf *.tar.gz
RUN cd /root/clustalw/clustalw* && cp * /usr/local/bin

# MUSCLE
ENV muscle_ftp http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux32.tar.gz
RUN mkdir /root/muscle && cd /root/muscle && wget ${muscle_ftp} && tar -zxvf * && rm -rf *.tar.gz
RUN cd /root/muscle/ && cp muscle* /usr/local/bin/muscle

# Python deps
RUN pip3 install natsort 
RUN pip3 install PyYAML
RUN apt-get update
RUN pip3 install matplotlib_venn
RUN apt-get update
RUN apt-get -y install python3-tk
RUN pip3 install rpy2
RUN pip3 install datrie
RUN pip3 install scikit-posthocs
RUN R -e "install.packages('FSA') "

