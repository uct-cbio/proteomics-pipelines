FROM cbio/bio:latest

ENV http_proxy ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

RUN apt-get update
RUN mkdir /root/fasta
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# Leave to end
RUN echo "export PATH=$PATH::/root/bin/python:/root/bin/R:/root/bin/bash" >> ~/.bashrc && source ~/.bashrc
RUN echo "export PYTHONPATH=$PYTHONPATH::/root/lib" >> ~/.bashrc && source ~/.bashrc

COPY lib /root/lib
COPY bin /root/bin

COPY docker/sixframe/config.yaml /root/config.yaml

CMD export PATH=$PATH::/root/bin/python:/root/bin/R:/root/bin/bash && export PYTHONPATH=$PYTHONPATH::/root/lib/ && sixframe.sh /root/fasta /root/config.yaml

