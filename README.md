EESARDocker
===========
Source files for the Docker image mort/eesardocker
Docker image with Jupyter widget interface for 
Sequential SAR change detection on Google Earth Engine

Pull and/or run the container with

    docker run -d -p 8888:8888 --name=eesar mort/eesardocker  
    
or from a raspberry pi with

	docker run -d -p 8888:8888 --name=eesar mort/rpi-eesardocker  
	
If you wish to link a local host working directory to the container, include the flag

	-v <path/to/mylocalfolder>:/home/mylocalfolder	

Point your browser to http://localhost:8888 to see the Jupyter notebook home page. 
Open a terminal window in the container from the notebook menu and enter

	earthengine authenticate
	
and follow the authentication instructions. This only has to be done once.	
 
Open the Notebook 

    interface.ipynb 

Stop the container with

    docker stop eesar 
     
Re-start with

    docker start eesar    