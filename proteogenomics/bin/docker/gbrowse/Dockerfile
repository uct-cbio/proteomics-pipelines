# InterProScan
# VERSION 0.1
# Tracking 5.21-60.0

FROM debian:jessie

# make sure the package repository is up to date
RUN DEBIAN_FRONTEND=noninteractive apt-get -qq update

RUN apt-get -y install gbrowse gbrowse-calign gbrowse-data libbio-samtools-perl apache2

#RUN apt-get -y install mysql-server

#RUN apt-get -y remove gbrowse

#RUN perl -MCPAN -e 'install Bio::Graphics::Browser2'

CMD /etc/init.d/apache2 start






