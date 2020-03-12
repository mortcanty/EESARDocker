#!/usr/bin/env python
#******************************************************************************
#  Name:     fidelity.py
#  Purpose:  
#    calculate the mean square bias between two polSAR images.
#
#  Usage:        
#    python fidelity.py [OPTIONS] filename1 filename2 
#
#  Copyright (c) 2018 Mort Canty

import numpy as np
import sys, getopt
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly
  
    
def main(): 
    usage = '''
    Usage:
------------------------------------------------

Calculate the mean square difference between two polSAR images   

python %s [OPTIONS] filename1 filename2
      
Options:

   -h          this help
   -p  <list>    band positions e.g. -p [1,2,3,4,5,7]
   -d <list>   spatial subset
'''      
   
    options, args = getopt.getopt(sys.argv[1:],'hd:p:') 
    dims = None 
    pos = None
    for option, value in options:
        if option == '-h':
            print(usage)
            return    
        elif option == '-d':
            dims = eval(value)   
        elif option == '-p':
            pos = eval(value)            
    if len(args)==2:
        fn1 = args[0] 
        fn2 = args[1]  
    else:
        print('Incorrect number of arguments')
        print(usage)
        sys.exit(1)        
    print( 'Comparing %s with %s'%(fn1,fn2))     
    gdal.AllRegister()
    inDataset = gdal.Open(fn1,GA_ReadOnly)  
    cols = inDataset.RasterXSize
    rows = inDataset.RasterYSize   
    bands = inDataset.RasterCount  
    if dims:
        x0,y0,cols,rows = dims
    else:
        x0 = 0
        y0 = 0    
    if pos is not None:
        bands = len(pos)
    else:
        pos = range(1,bands+1)            
    g1 = np.zeros((cols*rows,bands))    
    k = 0  
    for b in pos:      
        g1[:,k] = np.nan_to_num(inDataset.GetRasterBand(b).ReadAsArray(x0,y0,cols,rows).astype(float).ravel())
        k += 1       
    inDataset = gdal.Open(fn2,GA_ReadOnly)   
    g2 = np.zeros((cols*rows,bands))  
    k = 0 
    for b in pos:      
        g2[:,k] = np.nan_to_num(inDataset.GetRasterBand(b).ReadAsArray(x0,y0,cols,rows).astype(float).ravel()) 
        k += 1 
    print('minus log mean bias: %f' %-np.log(np.abs((np.mean(g2)-np.mean(g1))/np.mean(g1))))
    inDataset = None
   
if __name__ == '__main__':
    main()    