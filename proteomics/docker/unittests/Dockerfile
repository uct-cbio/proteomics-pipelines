FROM cbio/iprscan:latest

ENV http_proxy ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}
ENV PYTHONPATH $PYTHONPATH::/root/cbio-pipelines/proteomics/lib
ENV PATH $PATH::/root/cbio-pipelines/proteomics/bin/R:/root/cbio-pipelines/proteomics/bin/bash 
ENV python2ve /root/ve27
RUN echo "ready for unittests"
