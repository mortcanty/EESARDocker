#!/usr/bin/env python
#******************************************************************************
#  Name:     rt.py
#  Purpose:  
#    calculate the ratio of two polSAR images.
#
#  Usage:        
#    python ra.py [OPTIONS] filename1 filename2 
#
#  Copyright (c) 2020 Mort Canty

import numpy as np
import sys, os, getopt
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly, GDT_Float32
  
def main(): 
    usage = '''
    Usage:
------------------------------------------------

Calculate the ratio of two polSAR images   

python %s [OPTIONS] filename1 filename2
      
Options:

   -h          this help
   -p  <list>  band positions e.g. -p [1,2,3,4,5,7]
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
        
    path = os.path.dirname(fn1)    
    basename1 = os.path.basename(fn1)
    basename2 = os.path.basename(fn2)
    root1, ext = os.path.splitext(basename1)
    root2, ext = os.path.splitext(basename2)
    outfile = path + '/(' + root1 + '-over-' + root2 + ')' + ext        
    
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
        g1[:,k] = np.nan_to_num(inDataset.GetRasterBand(b).ReadAsArray(x0,y0,cols,rows).ravel())
        k += 1       
    inDataset = gdal.Open(fn2,GA_ReadOnly)   
    g2 = np.zeros((cols*rows,bands))  
    k = 0 
    for b in pos:      
        g2[:,k] = np.nan_to_num(inDataset.GetRasterBand(b).ReadAsArray(x0,y0,cols,rows).ravel()) 
 #       g2[:,k] = np.nan_to_num(inDataset.GetRasterBand(b).ReadAsArray(0,0,cols,rows).ravel()) 
        k += 1  
    
    g2 = np.where(g2==0,0.0001,g2)
    ratio = g1/g2


    driver = inDataset.GetDriver()   
    outDataset = driver.Create(outfile,cols,rows,bands,GDT_Float32)
    projection = inDataset.GetProjection()
    geotransform = inDataset.GetGeoTransform()
    if geotransform is not None:
        gt = list(geotransform)
        gt[0] = gt[0] + x0*gt[1]
        gt[3] = gt[3] + y0*gt[5]
        outDataset.SetGeoTransform(tuple(gt))
    if projection is not None:
        outDataset.SetProjection(projection)   
    for k in range(bands):       
        outBand = outDataset.GetRasterBand(k+1)
        outBand.WriteArray(np.reshape(ratio[:,k],(rows,cols)),0,0) 
        outBand.FlushCache() 
    outDataset = None  
    print( 'Ratio image written to: %s'%outfile )     
    
if __name__ == '__main__':
    main()    