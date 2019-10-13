# Dockerfile with Jupyter widget interface for 
# Sequential SAR on Earth Engine
 
FROM     debian:stretch

MAINTAINER Mort Canty "mort.canty@gmail.com"

ENV     REFRESHED_AT 2019-04-16
 
RUN apt-get update && apt-get install -y \
    python \
    build-essential \
    pandoc \
    python-dev \
    python-pygments \
    libssl-dev \
    python-pip  
    
RUN     pip install --upgrade pip   
RUN     apt-get install libffi-dev
COPY    requirements.txt /home 
RUN     pip install -r /home/requirements.txt       

EXPOSE 8888

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
ADD     scripts /home/scripts

# ipython notebook startup script
COPY    notebook.sh /
RUN     chmod u+x /notebook.sh

WORKDIR /home  
CMD     ["/notebook.sh"]
