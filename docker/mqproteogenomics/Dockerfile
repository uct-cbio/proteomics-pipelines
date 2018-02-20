FROM cbio/iprscan:latest

ENV http_proxy ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

RUN pip3 install matplotlib
RUN pip3 install rpy2
RUN R -e "install.packages('epitools', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('plotrix', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('sciplot', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('dunn.test', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('FSA',,'http://www.rforge.net/', dependencies=TRUE) "

RUN pip3 install urllib3
RUN pip3 install datrie

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# Leave to end
ENV PATH $PATH::/root/bin/python:/root/bin/R:/root/bin/bash
ENV PYTHONPATH $PYTHONPATH::/root/lib

COPY lib /root/lib
COPY bin /root/bin

#COPY docker/mqproteogenomics/config.yaml /root/config.yaml
CMD mqproteogenomics.sh /root/config.yaml /root/output 'mtu' /root/ve27 
#CMD export PATH=$PATH::/root/bin/python:/root/bin/R:/root/bin/bash && export PYTHONPATH=$PYTHONPATH::/root/lib/ && mqproteogenomics.sh /root/config.yaml

