FROM cbio/bio:latest

ARG http_proxy
ARG https_proxy
ARG ftp_proxy

ENV http_proxy  ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}


RUN pip3 install mygene

ENV DOCKER_JAVA_ARGS -XX:+UseParallelGC -XX:ParallelGCThreads= -XX:+AggressiveOpts -XX:+UseFastAccessorMethods -Xms128M -Xmx2048M
#COPY interproscan.properties /interproscan/interproscan.properties

RUN pip3 install tzlocal rpy2
ENV PYTHONPATH $PYTHONPATH::/home/bio/proteomics-pipelines/lib
ENV PATH $PATH::/home/bio/proteomics-pipelines/bin/R:/home/bio/proteomics-pipelines/bin/bash
ENV python2ve /home/bio/ve27
RUN virtualenv -p python2.7 /home/bio/ve27
