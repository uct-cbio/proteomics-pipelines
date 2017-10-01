FROM python:3.6-alpine

ENV http_proxy  ${http_proxy}

ENV https_proxy ${https_proxy}

RUN apk add --no-cache g++ && \
    ln -s /usr/include/locale.h /usr/include/xlocale.h && \
    pip install cython==0.25.2 numpy==1.12.0 && \
    pip install pandas==0.19.2 && \
    pip3 install elasticsearch 

# Bind a folder containing csv files to this folder (all files will be uploaded recursively
RUN mkdir /root/data

COPY bin/python/csv2es.py /root/csv2es.py

CMD python /root/csv2es.py ${ES_HOST} ${ES_PORT} ${ES_ALIAS} /root/data
