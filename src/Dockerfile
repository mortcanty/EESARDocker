# Dockerfile with Jupyter widget interface for 
# Sequential SAR on Earth Engine
 
FROM     debian:stretch

MAINTAINER Mort Canty "mort.canty@gmail.com"

ENV     REFRESHED_AT 2021-03-21
 
RUN apt-get update && apt-get install -y \
    python3 \
    build-essential \
    pandoc \
    python3-dev \
    python3-pygments \
    libssl-dev \
    libffi-dev \
    libpq-dev \
    gdal-bin \
    libgdal-dev \
    python3-gdal \
    python3-pip   
    
ENV     REFRESHED_AT 2022-06-24    
    
RUN     pip3 install --upgrade pip   

COPY    requirements.txt /home 
RUN     pip3 install -r /home/requirements.txt


# jupyter notebook with leaflet
RUN     pip3 install --upgrade pip
RUN     pip3 install numpy
RUN     pip3 install scipy
RUN     pip3 install jupyter
RUN     pip3 install ipyleaflet
RUN     jupyter nbextension enable --py --sys-prefix ipyleaflet



RUN     pip3 install -U cython

RUN     pip3 install matplotlib==3.0.3

EXPOSE 8888

# setup for earthengine
RUN     pip3 install pyasn1 --upgrade
RUN     pip3 install --upgrade setuptools && \
        pip3 install google-api-python-client && \
        pip3 install --upgrade oauth2client && \
        pip3 install pyCrypto && \
        apt-get install -y libssl-dev
        
RUN     pip3 install earthengine-api      
RUN     pip3 install earthengine-api --upgrade

#RUN     pip3 install --upgrade --no-cache-dir setuptools==57.0.0

#RUN     pip3 install --global-option=build_ext --global-option="-I/usr/include/gdal" GDAL==2.1.3


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

ENV     REFRESHED_AT 2022-06-23

COPY    interface.ipynb /home/interface.ipynb
COPY    dynearthnet.ipynb /home/dynearthnet.ipynb

ADD     scripts /home/scripts

# ipython notebook startup script
COPY    notebook.sh /
RUN     chmod u+x /notebook.sh

WORKDIR /home  
CMD     ["/notebook.sh"]
