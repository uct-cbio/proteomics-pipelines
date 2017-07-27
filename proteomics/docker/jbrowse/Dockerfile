FROM jbrowse/gmod-jbrowse:latest

ENV http_proxy  ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

ENV PATH $PATH::/jbrowse/bin
RUN apt-get update
RUN apt-get -y install nano vim
