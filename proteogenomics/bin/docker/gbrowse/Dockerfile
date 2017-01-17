# InterProScan
# VERSION 0.1
# Tracking 5.21-60.0

FROM ubuntu:latest

ARG http_proxy
ARG https_proxy

ENV http_proxy  ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

# make sure the package reposito

RUN apt-get update
RUN apt-get -y install wget make  
RUN apt-get -y install gbrowse gbrowse-calign gbrowse-data libbio-samtools-perl apache2
RUN apt-get -y remove gbrowse gbrowse-calign gbrowse-data
RUN cd /root && wget http://hgdownload.cse.ucsc.edu/admin/jksrc.zip && unzip jksrc.zip
RUN apt-get install -y libpng-dev 
RUN apt-get install -y gcc 
RUN apt-get install -y libssl-dev
RUN apt-get install -y wget subversion git-core libc6-dev build-essential libgd2-noxpm-dev
RUN apt-get install -y apache2 libapache2-mod-fcgid libapache2-mod-perl2 sqlite3 libdbd-sqlite3-perl
RUN apt-get install -y libmodule-build-perl libio-string-perl \
 libcapture-tiny-perl  libcgi-session-perl libwww-perl libstatistics-descriptive-perl \
 libjson-perl libjson-any-perl libsvg-perl perl-doc libtemplate-perl libgd-gd2-perl libgd-svg-perl \
 libdata-stag-perl libterm-readkey-perl libhttp-daemon-perl libfcgi-perl
RUN apt-get -y install postgresql libdbd-pg-perl
RUN apt-get -y install libdigest-sha-perl libssl-dev libmath-bigint-gmp-perl \
  libnet-openid-consumer-perl libnet-smtp-ssl-perl libauthen-sasl-perl libcrypt-ssleay-perl

RUN git clone https://github.com/bioperl/bioperl-live.git
RUN cd bioperl-live  && git checkout bioperl-release-1-6-901 && yes | perl Build.PL && yes | ./Build installdeps && ./Build install

ENV MACHTYPE x86_64 
RUN cd /root/kent/src/lib && make
ENV KENT_SRC /root/kent/src
RUN apt-get -y install git
RUN cd /root && git clone https://github.com/GMOD/GBrowse
RUN cd  /root/GBrowse  && perl Build.PL
RUN cd /root/GBrowse && yes | ./Build installdeps
RUN cd /root/GBrowse && ./Build test
RUN cd /root/GBrowse && ./Build apache_conf
RUN echo "ServerName localhost" | tee /etc/apache2/conf-available/fqdn.conf
RUN a2enconf fqdn
RUN /etc/init.d/apache2 restart
RUN apt-get -y install sudo

ENV http_proxy '' 
ENV https_proxy ''
ENV ftp_proxy ''

RUN adduser user 
RUN adduser user user
RUN echo 'user:user' | chpasswd
RUN cd /root/GBrowse && su user && ./Build install

CMD sudo service apache2 start && ln -sf /proc/self/fd/1 /var/log/apache2/access.log && \
    ln -sf /proc/self/fd/1 /var/log/apache2/error.log && tail -F /proc/self/fd/1
