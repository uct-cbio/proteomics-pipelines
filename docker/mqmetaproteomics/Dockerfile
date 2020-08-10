#FROM cbio/iprscan:latest

FROM thyscbio/bio:v1.2

RUN pip3 install xmltramp2

ENV http_proxy ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# Leave to end
ENV PATH $PATH::/home/bio/bin/python:/home/bio/bin/R:/home/bio/bin/bash:/home/bio/bin/python
ENV PYTHONPATH $PYTHONPATH::/home/bio/lib

COPY lib /home/bio/lib
COPY bin /home/bio/bin
