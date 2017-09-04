FROM cbio/bio:latest

ENV http_proxy ${http_proxy}
ENV https_proxy ${https_proxy}
ENV ftp_proxy ${ftp_proxy}

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# Leave to end
ENV PATH $PATH::/root/bin/python:/root/bin/R:/root/bin/bash:/root/bin/python
ENV PYTHONPATH $PYTHONPATH::/root/lib

COPY lib /root/lib
COPY bin /root/bin
CMD mqmetaproteomics.sh /root/txt/config.yaml
