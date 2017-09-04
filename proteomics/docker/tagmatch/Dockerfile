FROM cbio/bio:latest

ENV http_proxy ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

RUN apt-get update
RUN apt-get install -y wget tar
RUN wget http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.2.20/SearchGUI-3.2.20-mac_and_linux.tar.gz
RUN tar -zxvf *.tar.gz && rm -rf *.tar.gz
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz
RUN tar -zxvf *.tar.gz && rm -rf *.tar.gz

# Leave to end
ENV PATH $PATH::/root/bin/python:/root/bin/R:/root/bin/bash
ENV PYTHONPATH $PYTHONPATH::/root/lib

COPY lib /root/lib
COPY bin /root/bin

