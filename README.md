EESARDocker
===========
Source files for the Docker image mort/eesardocker
Docker image with Jupyter widget interface for 
Sequential SAR on Earth Engine

Pull and/or run the container with

    sudo docker run -d -p 443:8888 --name=eesar mort/eesardocker 
    
or, on the Rasberry Pi,     

    sudo docker run -d -p 443:8888 --name=eesar mort/rpi-eesardocker

Point your browser to http://localhost:443 to see the Jupyter notebook home page. 

Open the Notebook 

    interface.ipynb 

Stop the container with

    sudo docker stop eesar 
     
Re-start with

    sudo docker start eesar    