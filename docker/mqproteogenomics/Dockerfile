FROM thyscbio/bio:v1.3

ENV http_proxy ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}


#COPY docker/mqproteogenomics/config.yaml /root/config.yaml
CMD mqproteogenomics.sh /root/config.yaml /root/output 'mtu' /root/ve27 
#CMD export PATH=$PATH::/root/bin/python:/root/bin/R:/root/bin/bash && export PYTHONPATH=$PYTHONPATH::/root/lib/ && mqproteogenomics.sh /root/config.yaml

