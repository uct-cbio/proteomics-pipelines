sudo docker build --build-arg http_proxy=${http_proxy} --build-arg https_proxy=$https_proxy -t cbio/jbrowse:latest .
sudo docker run -it --rm -p 8088:80 -v $1:/root/data --name cbio_jbrowse cbio/jbrowse:latest 
