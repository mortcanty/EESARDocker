# Dockerfile with Jupyter widget interface for 
# Sequential SAR on Earth Engine
 
FROM     debian:jessie  

# arm (raspberry pi) architecture       
#FROM     resin/rpi-raspbian:jessie 
           
MAINTAINER Mort Canty "mort.canty@gmail.com"
ENV     REFRESHED_AT 2018-05-17

RUN apt-get update && apt-get install -y \
    python \
    build-essential \
    make \
    gcc \
    pandoc \
    python-dev \
    python-pygments \
    python-pip  

# ipython notebook
RUN     pip install --upgrade pip
RUN     pip install jupyter
RUN     pip install ipyleaflet
RUN     jupyter nbextension enable --py --sys-prefix ipyleaflet


EXPOSE 8888

# setup for earthengine
RUN     pip install pyasn1 --upgrade
RUN     pip install --upgrade setuptools && \
        pip install google-api-python-client && \
        pip install --upgrade oauth2client && \
        pip install pyCrypto && \
        apt-get install -y libssl-dev
RUN     pip install earthengine-api

# install auxil
COPY    dist/auxil-1.0.tar.gz /home/auxil-1.0.tar.gz
WORKDIR /home
RUN     tar -xzvf auxil-1.0.tar.gz
WORKDIR /home/auxil-1.0
RUN     python setup.py install  
WORKDIR /home
RUN     rm -rf auxil-1.0
RUN     rm auxil-1.0.tar.gz

COPY    interface.ipynb /home/interface.ipynb

# ipython notebook startup script
COPY    notebook.sh /
RUN     chmod u+x /notebook.sh

WORKDIR /home  
CMD     ["/notebook.sh"]
