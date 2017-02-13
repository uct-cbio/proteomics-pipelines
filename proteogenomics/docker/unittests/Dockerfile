FROM cbio/bio:latest

ENV http_proxy ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

RUN apt-get update
ENV PYTHONPATH $PYTHONPATH::/root/cbio-pipelines/proteogenomics/lib

RUN echo "ready for unittests"
