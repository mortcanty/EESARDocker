# Dockerfile with Jupyter widget interface for 
# Sequential SAR on Earth Engine
 
FROM     debian:latest

MAINTAINER Mort Canty "mort.canty@gmail.com"

ENV     REFRESHED_AT 2020-01-16
 
RUN apt-get update && apt-get install -y \
    python3 \
    build-essential \
    pandoc \
    python3-dev \
    python3-pygments \
    libssl-dev \
    libffi-dev \
    gdal-bin \
    libgdal-dev \
    python3-pip  
    
RUN     pip3 install --upgrade pip   

EXPOSE 8888

COPY    requirements.txt /home 
RUN     pip3 install -r /home/requirements.txt
RUN     pip3 install --global-option=build_ext --global-option="-I/usr/include/gdal" GDAL==2.1.3

RUN     jupyter nbextension enable --py --sys-prefix ipyleaflet

# setup the prov_means library
COPY    prov_means.c /home/prov_means.c
WORKDIR /home
RUN     gcc -shared -Wall -g -o libprov_means.so -fPIC prov_means.c
RUN     mv libprov_means.so /usr/lib/
RUN     rm prov_means.c

# install auxil
COPY    dist/auxil-1.0.tar.gz /home/auxil-1.0.tar.gz
WORKDIR /home
RUN     tar -xzvf auxil-1.0.tar.gz
WORKDIR /home/auxil-1.0
RUN     python3 setup.py install  
WORKDIR /home
RUN     rm -rf auxil-1.0
RUN     rm auxil-1.0.tar.gz

COPY    interface.ipynb /home/interface.ipynb
# COPY    notebooks/s1class.ipynb /home/s1class.ipynb

ADD     scripts /home/scripts

# ipython notebook startup script
COPY    notebook.sh /
RUN     chmod u+x /notebook.sh

WORKDIR /home  
CMD     ["/notebook.sh"]
