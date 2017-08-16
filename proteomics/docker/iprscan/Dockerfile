FROM cbio/bio:latest

ARG http_proxy
ARG https_proxy
ARG ftp_proxy

ENV http_proxy  ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

ENV ipr_ftp "ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.21-60.0/interproscan-5.21-60.0-64-bit.tar.gz"
ENV panther_ftp "ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-10.0.tar.gz"

# make sure the package repository is up to date
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install software-properties-common

# Install Oracle Java 8
RUN \
  echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
  add-apt-repository -y ppa:webupd8team/java && \
  apt-get update && \
  apt-get install -y oracle-java8-installer && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /var/cache/oracle-jdk8-installer
ENV JAVA_HOME /usr/lib/jvm/java-8-oracle

# 1. Obtaining the core InterProScan software
RUN \
  wget ${ipr_ftp}  && \
  wget ${ipr_ftp}.md5 && \
  md5sum -c $( basename ${ipr_ftp}.md5) && \
  tar -pxvzf $(basename ${ipr_ftp}) && \
  rm $(basename ${ipr_ftp}) && \
  rm $( basename ${ipr_ftp}.md5)

# 2. Installing Panther Models
RUN \
  cd /interproscan-*/data && \
  wget ${panther_ftp} && wget ${panther_ftp}.md5

RUN \
  cd /interproscan-*/data && \
  md5sum -c $( basename ${panther_ftp}.md5) 

RUN \
  cd /interproscan-*/data && \
  tar -pxvzf $( basename ${panther_ftp})

RUN mv /interproscan-* /interproscan
ENV PATH ${PATH}::/interproscan

#VOLUME ["/interproscan-5.21-60.0/data/"]
#COPY ./interproscan.properties /interproscan-5.21-60.0/interproscan.properties
#COPY ./interproscan.sh /interproscan-5.21-60.0/interproscan.sh

ENV DOCKER_JAVA_ARGS -XX:+UseParallelGC -XX:ParallelGCThreads= -XX:+AggressiveOpts -XX:+UseFastAccessorMethods -Xms128M -Xmx2048M

