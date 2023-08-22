#FROM thyscbio/bio:v2.10
FROM thyscbio/bio:v3.0.17

# samtools etc
#RUN apt-get update && apt-get install -y libssl-dev openssl
#RUN conda update conda && conda install -c bioconda bcftools minimap2
#RUN conda config --add channels bioconda && conda config --add channels conda-forge 
#RUN conda install samtools 
#RUN conda install tabix 
#RUN conda install jbrowse2
#RUN conda install emboss
#RUN yes | conda update --all
#RUN conda init
#RUN conda activate 
# install jbrowse2
#RUN apt-get update && apt-get install -y \
#    software-properties-common build-essential zlib1g-dev  genometools
#RUN npm install npm@latest -g && \
#    npm install n -g && \
#    n latest
#RUN npm install -g @jbrowse/cli
#RUN npm install -g serve
#RUN jbrowse create jbrowse_folder
COPY lib /home/bio/lib 
COPY bin /home/bio/bin

#COPY docker/mqproteogenomics/config.yaml /root/config.yaml
CMD mqproteogenomics.sh /root/config.yaml /root/output 'mtu' /root/ve27 
#CMD export PATH=$PATH::/root/bin/python:/root/bin/R:/root/bin/bash && export PYTHONPATH=$PYTHONPATH::/root/lib/ && mqproteogenomics.sh /root/config.yaml

