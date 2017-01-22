
FROM ubuntu:latest

# updated 20 Jan 2017

ARG http_proxy
ARG https_proxy

ENV http_proxy  ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

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
RUN apt-get -y install postgresql libdbd-pg-perl sudo
RUN apt-get -y install libdigest-sha-perl libssl-dev libmath-bigint-gmp-perl \
  libnet-openid-consumer-perl libnet-smtp-ssl-perl libauthen-sasl-perl libcrypt-ssleay-perl

RUN git clone https://github.com/bioperl/bioperl-live.git
RUN cd bioperl-live  && git checkout bioperl-release-1-6-901 && yes | perl Build.PL && yes | ./Build installdeps && ./Build install


#MySql
RUN echo "mysql-server mysql-server/root_password password gbrowse2" | sudo debconf-set-selections
RUN echo "mysql-server mysql-server/root_password_again password gbrowse2" | sudo debconf-set-selections
RUN apt-get update && apt-get -y install mysql-client
RUN apt-get update && apt-get -y install mysql-server
RUN apt-get update && apt-get -y install libmysqlclient-dev

#Kent 
ENV MACHTYPE x86_64 
RUN cd /root/kent/src/lib && make
RUN cd /root/kent/src/jkOwnLib && make
RUN cd /root/kent/src/hg/lib && make
RUN cd /root/kent/src/htslib && make
RUN mkdir /root/bin && mkdir /root/bin/x86_64
RUN cd /root/kent/src/hg/ratStuff/mafFilter && make
ENV KENT_SRC /root/kent/src

#Gbrowse
RUN cd /root && git clone https://github.com/GMOD/GBrowse
RUN cd  /root/GBrowse  && perl Build.PL
RUN cd /root/GBrowse && yes | ./Build installdeps
RUN cd /root/GBrowse && ./Build test
RUN cd /root/GBrowse && ./Build apache_conf
RUN echo "ServerName localhost" | tee /etc/apache2/conf-available/fqdn.conf
RUN a2enconf fqdn
RUN /etc/init.d/apache2 restart
RUN apt-get -y install sudo

#Config files and text editors
RUN apt-get -y install vim nano
RUN cd /root && git clone https://thys_potgieter@bitbucket.org/thys_potgieter/configs.git
RUN cp /root/configs/.vimrc /root/.vimrc
RUN cp /root/configs/.tmux.conf /root/.tmux.conf
RUN apt-get install -y libdbd-sqlite3-perl cmake groff info python3 python3-pip python3-dev
RUN pip3 install numpy scipy pandas 
RUN pip3 install BioPython 

RUN adduser user 
RUN adduser user user
RUN echo 'user:user' | chpasswd
RUN cd /root/GBrowse && su user && ./Build install


# install mugsy
ENV MUGSY_INSTALL "/root/mugsy_x86-64-v1r2.3"
# below is just if sourceforge blocks you 
COPY mugsy_x86-64-v1r2.3.tgz /root/download
#RUN cd /root && wget https://sourceforge.net/projects/mugsy/files/mugsy_x86-64-v1r2.3.tgz/download
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN cd /root && tar xvzf download && cd mugsy* && echo "export PATH=$PATH::"`pwd` >> ~/.bashrc && echo "source "`pwd`"/mugsyenv.sh" >> ~/.bashrc

#MafFilter
RUN cd /root && wget http://biopp.univ-montp2.fr/repos/exe/lin64/maffilter/maffilter-x86_64-bin-static-1.1.2-1.tar.gz
RUN cd /root && tar -zxvf maffilter*.gz
RUN cd /root && rm -rf *.gz 
RUN echo "export PATH=$PATH:/root/MafFilter:/root/mugsy_x86-64-v1r2.3" >> ~/.bashrc


# Persistent volumes
VOLUME  ["/etc/gbrowse2", "/var/lib/gbrowse2/databases", "/etc/mysql/", "/var/lib/mysql"]

RUN echo "export PATH=$PATH:/root/MafFilter:/root/mugsy_x86-64-v1r2.3:/root/bin/x86_64" >> ~/.bashrc
COPY SyntenyIO.pm /usr/local/lib/x86_64-linux-gnu/perl/5.22.1/Legacy/DB/SyntenyIO.pm
COPY syntenydata.py /root/syntenydata.py
COPY mugsyenv.sh $MUGSY_INSTALL

RUN cd /var/www && ln -s ../tmp
RUN mkdir /var/www/tmp
RUN mkdir /var/www/tmp/gbrowse2
RUN chmod 777 /var/www/tmp/gbrowse2

# Remove proxy
ENV http_proxy '' 
ENV https_proxy ''
ENV ftp_proxy ''
CMD sudo service apache2 start && sudo service mysql start && ln -sf /proc/self/fd/1 /var/log/apache2/access.log && \
    ln -sf /proc/self/fd/1 /var/log/apache2/error.log && tail -F /proc/self/fd/1
