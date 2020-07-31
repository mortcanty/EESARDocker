#!/usr/bin/env python
#******************************************************************************
#  Name:     atsfthreshold.py
#  Purpose:  
#    calculate threshold ATSF with frequency image and replace with refined Lee
#
#  Usage:        
#    python atsfthreshold.py [OPTIONS] atsfFile leefile avimgLogFile
#
#  Copyright (c) 2018 Mort Canty

import numpy as np
import sys, getopt, os
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly,GDT_Float32
  
    
def main(): 
    usage = '''
    Usage:
------------------------------------------------

Threshold ATSF with frequency image and replace with refined Lee filter  

python %s [OPTIONS] atsffile leeFile avimglogfile 
      
Options:

   -h          this help
   -d <list>   spatial subset
   -t <int>    threshold (default maxfreq/4)
'''      
   
    options, args = getopt.getopt(sys.argv[1:],'hd:t:') 
    dims = None 
    thresh = None
    for option, value in options:
        if option == '-h':
            print(usage)
            return    
        elif option == '-d':
            dims = eval(value)   
        elif option == '-t':
            thresh = eval(value)  
    if len(args)==3:
        fn1 = args[0] 
        fn2 = args[1]  
        fn3 = args[2]
    else:
        print('Incorrect number of arguments')
        print(usage)
        sys.exit(1)        
        
    path = os.path.dirname(fn1)    
    basename = os.path.basename(fn1)
    root, ext = os.path.splitext(basename)
    outfn = path + '/' + root + '_thresh' + ext
        
    print( 'Replacing %s with %s under threshold'%(fn1,fn2))     
    gdal.AllRegister()
    inDataset1 = gdal.Open(fn1,GA_ReadOnly)  
    cols = inDataset1.RasterXSize
    rows = inDataset1.RasterYSize   
    bands = inDataset1.RasterCount  
    inDataset2 = gdal.Open(fn2,GA_ReadOnly)  
    cols2 = inDataset2.RasterXSize
    rows2 = inDataset2.RasterYSize   
    cols = min(cols,cols2) 
    rows = min(rows,rows2)
    if dims:
        x0,y0,cols,rows = dims
    else:
        x0 = 0
        y0 = 0    
    g1 = np.zeros((cols*rows,bands))      
    for k in range(bands):      
        g1[:,k] = np.nan_to_num(inDataset1.GetRasterBand(k+1).ReadAsArray(x0,y0,cols,rows).astype(float).ravel())       
    g2 = np.zeros((cols*rows,bands))   
    for k in range(bands):      
        g2[:,k] = np.nan_to_num(inDataset2.GetRasterBand(k+1).ReadAsArray(x0,y0,cols,rows).astype(float).ravel())   
    inDataset = gdal.Open(fn3,GA_ReadOnly)   
    g3 = np.zeros((cols*rows,bands))   
    for k in range(bands):      
        g3[:,k] = np.nan_to_num(inDataset.GetRasterBand(1).ReadAsArray(x0,y0,cols,rows).astype(float).ravel())   
    if thresh == None:
        thresh = np.max(g3)/4   
    inDataset2 = None
    inDataset = None
    
    g1 = np.reshape(np.where(g3<thresh,g2,g1),(rows,cols,bands)) 
    
    driver = inDataset1.GetDriver() 
    outDataset = driver.Create(outfn,cols,rows,bands,GDT_Float32)
    
    geotransform = inDataset1.GetGeoTransform()
    if geotransform is not None:
        gt = list(geotransform)
        gt[0] = gt[0] + x0*gt[1]
        gt[3] = gt[3] + y0*gt[5]
        outDataset.SetGeoTransform(tuple(gt))
    projection = inDataset1.GetProjection()        
    if projection is not None:
        outDataset.SetProjection(projection) 
    for k in range(bands):    
        outBand = outDataset.GetRasterBand(k+1)
        outBand.WriteArray(g1[:,:,k],0,0) 
        outBand.FlushCache() 
    outDataset = None
    inDataset1 = None
    print('result written to: '+outfn) 
   
if __name__ == '__main__':
    main()    