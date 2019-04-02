FROM cbio/bio:latest

ARG http_proxy
ARG https_proxy
ARG ftp_proxy

ENV http_proxy  ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}
ENV ipr_ftp "ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.30-69.0/interproscan-5.30-69.0-64-bit.tar.gz"
ENV panther_ftp "ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-12.0.tar.gz"

# 1. Obtaining the core InterProScan software
#RUN \
#  wget ${ipr_ftp}  && \
#  wget ${ipr_ftp}.md5 && \
#  md5sum -c $( basename ${ipr_ftp}.md5) && \
#  tar -pxvzf $(basename ${ipr_ftp}) && \
#  rm $(basename ${ipr_ftp}) && \
#  rm $( basename ${ipr_ftp}.md5)

# 2. Installing Panther Models
#RUN \
#  cd /interproscan-*/data && \
#  wget ${panther_ftp} && wget ${panther_ftp}.md5

#RUN \
#  cd /interproscan-*/data && \
#  md5sum -c $( basename ${panther_ftp}.md5) 

#RUN \
#  cd /interproscan-*/data && \
#  tar -pxvzf $( basename ${panther_ftp})

#RUN mv /interproscan-* /interproscan
#RUN ls
#ENV PATH ${PATH}::/interproscan-5.26-65.0  
#RUN ln -s /interproscan-* /interproscan

#RUN cd /interproscan && wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz && tar xvf ncbi-blast-2.6.0+-x64-linux.tar.gz
#RUN cd /interproscan && cp bin/blast/ncbi-blast-2.6.0+/rpsblast bin/blast/ncbi-blast-2.6.0+/rpsblast.bak
#RUN cd /interproscan && cp ncbi-blast-2.6.0+/bin/rpsblast bin/blast/ncbi-blast-2.6.0+/rpsblast

#ENV PATH ${PATH}::/interproscan 
# make sure the package repository is up to date
#RUN apt-get update
#RUN DEBIAN_FRONTEND=noninteractive apt-get -y install software-properties-common

# Install Oracle Java 8
#RUN wget --header "Cookie: oraclelicense=accept-securebackup-cookie" http://download.oracle.com/otn-pub/java/jdk/8u162-b12/0da788060d494f5095bf8624735fa2f1/jdk-8u162-linux-x64.tar.gz

#RUN mkdir /opt/jdk
#RUN tar -zxf jdk-*x-x64.tar.gz -C /opt/jdk
#RUN ls /opt/jdk 
#RUN update-alternatives --install /usr/bin/java java /opt/jdk/jdk1.8.0_162/bin/java 100 
#RUN update-alternatives --install /usr/bin/javac javac /opt/jdk/jdk1.8.0_162/bin/javac 100 
#RUN update-alternatives --display java
#RUN update-alternatives --display javac
#RUN apt-get -y install default-jdk

#VOLUME ["/interproscan-5.21-60.0/data/"]
#COPY ./interproscan.properties /interproscan-5.21-60.0/interproscan.properties

#COPY ./interproscan.sh /interproscan-5.21-60.0/interproscan.sh
RUN pip3 install mygene

ENV DOCKER_JAVA_ARGS -XX:+UseParallelGC -XX:ParallelGCThreads= -XX:+AggressiveOpts -XX:+UseFastAccessorMethods -Xms128M -Xmx2048M
#COPY interproscan.properties /interproscan/interproscan.properties

RUN pip3 install tzlocal rpy2
ENV PYTHONPATH $PYTHONPATH::/root/proteomics-pipelines/lib
ENV PATH $PATH::/root/proteomics-pipelines/bin/R:/root/proteomics-pipelines/bin/bash
ENV python2ve /root/ve27

