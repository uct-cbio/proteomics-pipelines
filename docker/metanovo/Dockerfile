FROM thyscbio/bio:v2.2

ENV http_proxy ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

RUN apt-get update
RUN apt-get install -y wget tar bc parallel 


#RUN groupadd -r bio && useradd -r -g bio bio

#WORKDIR /home/bio

COPY bin/python/chunkfasta.py /home/bio/chunkfasta.py
# Get UniProt
#RUN wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

# SearchGUI
RUN wget http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.2.20/SearchGUI-3.2.20-mac_and_linux.tar.gz
#RUN wget http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/4.1.7/SearchGUI-4.1.7-mac_and_linux.tar.gz
RUN tar -zxvf *.tar.gz && rm -rf *.tar.gz
#RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz
#RUN tar -zxvf *.tar.gz && rm -rf *.tar.gz

# Compomics utilities
#RUN wget http://genesis.ugent.be/maven2/com/compomics/utilities/4.11.19/utilities-4.11.19.zip
#RUN unzip utilities-4.11.19.zip && rm -rf *.zip



RUN wget http://genesis.ugent.be/maven2/com/compomics/utilities/4.12.0/utilities-4.12.0.zip
#RUN wget http://genesis.ugent.be/maven2/com/compomics/utilities/5.0.38/utilities-5.0.38.zip
RUN unzip utilities-4.12.0.zip && rm -rf *.zip

# DeNovoGUI
RUN wget http://genesis.ugent.be/maven2/com/compomics/denovogui/DeNovoGUI/1.15.11/DeNovoGUI-1.15.11-mac_and_linux.tar.gz
RUN tar -zxvf *.tar.gz && rm -rf *.tar.gz

RUN apt-get update && apt-get -y install git cmake
RUN apt-get update && apt-get -y  install qtbase5-dev
RUN apt-get update && apt-get -y install qtdeclarative5-dev


ENV PATH="$PATH::/home/bio/bin/python:/home/bio/bin/R:/home/bio/bin/bash"

ENV PYTHONPATH="$PYTHONPATH::/home/bio/lib"
ENV SG_PATH="/home/bio/SearchGUI*"
ENV DG_PATH="/home/bio/DeNovoGUI*"
ENV CU_PATH="/home/bio/utilities*"


COPY lib /home/bio/lib
COPY bin /home/bio/bin
COPY docker/metanovo/config.sh /home/bio/config.sh
COPY docker/metanovo/run.sh /home/bio/run.sh
COPY docker/metanovo/sshloginfile /home/bio/sshloginfile

ENV PATH="/home/bio/bin/bash:/home/bio/bin/python:${PATH}"

