#!/usr/bin/env python
#******************************************************************************
#  Name:     envi2tif.py
#  Purpose:  
#    read s1 image bands in envi format, write as GeoTiff.
#
#  Usage:        
#    python envi2tig.py [OPTIONS] filename stub 
#
#  Copyright (c) 2020 Mort Canty

import numpy as np
import sys, getopt
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly, GDT_Float32
  
def main(): 
    usage = '''
    Usage:
------------------------------------------------

Read s1 image bands in envi format, write as GeoTiff.

python %s [OPTIONS]  filename stub
      
Options:

   -h          this help

'''      
   
    options, args = getopt.getopt(sys.argv[1:],'h') 
    for option, value in options:
        if option == '-h':
            print(usage)
            return    
    if len(args)==1:
        stub = args[0] 
    else:
        print('Incorrect number of arguments')
        print(usage)
        sys.exit(1)        
        
    
    gdal.AllRegister()
    
    C11 = stub+'_vvvv_ml_gec.envi'
    C12 = stub+'_vvvh_ml_gec.envi'
    C22 = stub+'_vhvh_ml_gec.envi'
    
    inDataset = gdal.Open(C11,GA_ReadOnly)  
    cols = inDataset.RasterXSize
    rows = inDataset.RasterYSize   
    
    G = np.zeros((cols*rows,4))
    
    G[:,0] = np.nan_to_num(inDataset.GetRasterBand(1).ReadAsArray(0,0,cols,rows).ravel())
    
    inDataset = gdal.Open(C22,GA_ReadOnly)
    G[:,3] = np.nan_to_num(inDataset.GetRasterBand(1).ReadAsArray(0,0,cols,rows).ravel())
    

    inDataset = gdal.Open(C12,GA_ReadOnly)
    D = np.nan_to_num(inDataset.GetRasterBand(1).ReadAsArray(0,0).ravel())
    
    G[:,1] = np.real(D)
    G[:,2] = np.imag(D)
    
    outfile = stub+'.tif'
    
    driver = gdal.GetDriverByName('GTiff')   
    outDataset = driver.Create(outfile,cols,rows,4,GDT_Float32)
    projection = inDataset.GetProjection()
    geotransform = inDataset.GetGeoTransform()
    if geotransform is not None:
        outDataset.SetGeoTransform(geotransform)
    if projection is not None:
        outDataset.SetProjection(projection)   
    for k in range(4):       
        outBand = outDataset.GetRasterBand(k+1)
        outBand.WriteArray(np.reshape(G[:,k],(rows,cols)),0,0) 
        outBand.FlushCache() 
    outDataset = None  
    print( 'GeoTiff image written to: %s'%outfile )     
    
if __name__ == '__main__':
    main()    