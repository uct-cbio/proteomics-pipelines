FROM cbio/bio:latest

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
RUN tar -zxvf *.tar.gz && rm -rf *.tar.gz
#RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz
#RUN tar -zxvf *.tar.gz && rm -rf *.tar.gz

# Compomics utilities
RUN wget http://genesis.ugent.be/maven2/com/compomics/utilities/4.11.19/utilities-4.11.19.zip
RUN unzip utilities-4.11.19.zip && rm -rf *.zip

# DeNovoGUI
RUN wget http://genesis.ugent.be/maven2/com/compomics/denovogui/DeNovoGUI/1.15.11/DeNovoGUI-1.15.11-mac_and_linux.tar.gz
RUN tar -zxvf *.tar.gz && rm -rf *.tar.gz

# PeptideShaker
RUN wget http://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/1.16.13/PeptideShaker-1.16.13.zip
RUN unzip PeptideShaker-*.zip && rm *.zip

# MSGFPlus

#RUN mkdir /msgfplus && cd /msgfplus && wget https://omics.pnl.gov/sites/default/files/MSGFPlus.zip && unzip MSGFPlus.zip
RUN apt-get update && apt-get -y install git cmake
RUN apt-get update && apt-get -y  install qtbase5-dev
RUN apt-get update && apt-get -y install qtdeclarative5-dev


#RUN mkdir /openms && cd /openms && git clone  https://github.com/OpenMS/contrib.git
#RUN cd /openms && mkdir contrib-build && cd contrib-build && cmake -DBUILD_TYPE=ALL -DNUMBER_OF_JOBS=4 ../contrib 
#RUN cd /openms && git clone https://github.com/OpenMS/OpenMS.git 
#RUN apt-get -y install build-essential cmake autoconf patch libtool automake
#RUN apt-get -y install qt4-default libqtwebkit-dev
#RUN apt-get -y install libeigen3-dev libwildmagic-dev libxerces-c-dev libboost-all-dev libsvn-dev libgsl-dev libbz2-dev
#RUN cd /openms && mkdir OpenMS-build && cd OpenMS-build && cmake -DOPENMS_CONTRIB_LIBS="/openms/contrib-build" -DBOOST_USE_STATIC=OFF ../OpenMS && make 

RUN wget http://www.proteoannotator.org/datasets/releases/ProteoAnnotator-1.7.86.zip && unzip ProteoAnnotator-*.zip


#ENV LD_LIBRARY_PATH /openms/OpenMS-build/lib:$LD_LIBRARY_PATH
ENV PATH="$PATH::/home/bio/bin/python:/home/bio/bin/R:/home/bio/bin/bash"

ENV PYTHONPATH="$PYTHONPATH::/home/bio/lib"
ENV SG_PATH="/home/bio/SearchGUI*"
ENV DG_PATH="/home/bio/DeNovoGUI*"
ENV CU_PATH="/home/bio/utilities*"
ENV TANDEM_DEFAULT_INPUT_PATH="/home/bio/default_input.xml"
ENV TANDEM_INPUT_STYLE_PATH="/home/bio/tandem-input-style.xsl"
ENV MZIDLIB_PATH="/home/bio/ProteoAnnotator-*/mzidlib*"

COPY docker/metanovo/tandem-input-style.xsl /home/bio/tandem-input-style.xsl
COPY docker/metanovo/default_input.xml /home/bio/default_input.xml

COPY lib /home/bio/lib
COPY bin /home/bio/bin
COPY docker/metanovo/config.sh /home/bio/config.sh
COPY docker/metanovo/run.sh /home/bio/run.sh
COPY docker/metanovo/sshloginfile /home/bio/sshloginfile

ENV PATH="/home/bio/bin/bash:/home/bio/bin/python:${PATH}"

# Install mono
#RUN apt-get -y install mono-complete

# Run denovo
#CMD /root/bin/bash/metanovo.sh  2>&1 | tee /root/output/log.txt 
