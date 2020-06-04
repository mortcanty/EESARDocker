#!/usr/bin/env python
#******************************************************************************
#  Name:     eemd_filter.py
#  Purpose:  Perform eemd spectral filtering on multi-temporal SAR intensity imagery        
#
#  Usage:             
#    python eemd_filetr.py [OPTIONS] filedir target
#
# MIT License
# 
# Copyright (c) 2020 Mort Canty


def getimg(fn,b,dims):
#  read band b of SAR polarimetric matrix file 
    from osgeo.gdalconst import GA_ReadOnly
    from osgeo import gdal
    import numpy as np
    import sys
    
    gdal.AllRegister()
    x0,y0,cols,rows = dims
    try:            
        inDataset = gdal.Open(fn,GA_ReadOnly)                     
        result = np.nan_to_num(inDataset.GetRasterBand(b+1).ReadAsArray(x0,y0,cols,rows).ravel())
        result = np.where(result<=0,10e-9,result)
        inDataset = None     
        return result  
    except Exception as e:
        print( 'Error: %s  -- Could not read file'%e )
        sys.exit(1)    
                       
def main():  
    import numpy as np
    import os, sys, time, getopt
    from os import listdir
    from os.path import isfile, join
    from osgeo import gdal 
    from osgeo.gdalconst import GA_ReadOnly, GDT_Float32
    from tempfile import NamedTemporaryFile
    from pyeemd import emd
    
    usage = '''
Usage:
------------------------------------------------

EEMD speckle filter for polarimetric SAR images

python %s [OPTIONS]  infiledir target 

Options:
  
  -h           this help 
  -d  <list>   spatial subset

infiledir:

  path to input files directory
  
target (<int>):

  index of target file to filter 
  
enl (<float>):

  equivalent number of looks

-------------------------------------------------'''%sys.argv[0]

    options,args = getopt.getopt(sys.argv[1:],'hd:')
    dims = None
    for option, value in options: 
        if option == '-h':
            print( usage )
            return 
        elif option == '-d':
            dims = eval(value)            
    if len(args)!=2:
        print('incorrect number of arguments')
        print( usage )
        sys.exit()
    pth = args[0]    
    fns1 = sorted([join(pth,f) for f in listdir(pth)])
    fns = []
    for f in fns1:
        if f.find('emd') == -1:
            fns.append(f)   
    k = len(fns)
    target = eval(args[1]) 
    gdal.AllRegister()        
#  first SAR image   
    try:            
        inDataset1 = gdal.Open(fns[0],GA_ReadOnly)                             
        cols = inDataset1.RasterXSize
        rows = inDataset1.RasterYSize    
        bands = inDataset1.RasterCount
    except Exception as e:
        print( 'Error: %s  -- Could not read file'%e)
        sys.exit(1)    
    if dims == None:
        dims = [0,0,cols,rows]
    x0,y0,cols,rows = dims            
    print( '===============================================' )
    print( '     EEMD Spectral filtering' )
    print( '===============================================' )   
    print( time.asctime() )  
    print( 'First image:  %s'%fns[0] )
    print( 'number of images: %i'%k )
    if bands==9:
        print( 'Quad polarization')
        pos = [0,5,8]
    elif bands==4:
        print( 'Dual polarizaton' )
        pos = [0,3]
    elif bands==3:
        print( 'Quad polarization, diagonal only' )
        pos = [0,1,2]
    elif bands==2:
        print( 'Dual polarization, diagonal only' )
        pos = [0,1]
    else:
        print( 'Intensity image' )
        pos = [0]
#  output file
    path = os.path.dirname(fns[0])    
    basename = os.path.basename(fns[target-1])
    print( 'Target image: %s'%basename )
    root, ext = os.path.splitext(basename)
    outfn = path + '/' + root + '_emd' + ext
#  create temporary, memory-mapped array of images (1st band only, in decibels)
    mm = NamedTemporaryFile()
    imarray = np.memmap(mm.name,dtype=np.float64,mode='w+',shape=(rows*cols,k))  
    outarr = np.zeros((len(pos),cols*rows))
    start = time.time()   
    ell = 0
    for p in pos:
        print( 'filtering band: %i'%p )       
        for j in range(k):
            arr = getimg(fns[j],p,dims)
            imarray[:,j] = 10*np.log10(arr)
#      run the ceemdan algorithm        
        for i in range(rows*cols):
            imfs = emd(imarray[i,:], S_number=4, num_siftings=50)
            imarray[i,:] = np.sum(imfs[2:,:],0)
#  restore linear scale to filtered target        
        outarr[ell,:] = 10**(imarray[:,target-1]/10.0)
        ell += 1
#  write to file system     
    driver = inDataset1.GetDriver() 
    outDataset = driver.Create(outfn,cols,rows,len(pos),GDT_Float32)
    geotransform = inDataset1.GetGeoTransform()
    if geotransform is not None:
        gt = list(geotransform)
        gt[0] = gt[0] + x0*gt[1]
        gt[3] = gt[3] + y0*gt[5]
    outDataset.SetGeoTransform(tuple(gt))
    projection = inDataset1.GetProjection()        
    if projection is not None:
        outDataset.SetProjection(projection)  
    ell = 0
    for p in pos:
        outBand = outDataset.GetRasterBand(ell+1)
        outBand.WriteArray(np.reshape(outarr[ell,:],(rows,cols)),0,0) 
        outBand.FlushCache()  
        ell += 1
    print('result written to: '+outfn) 
    print('elapsed time: '+str(time.time()-start))
    outDataset = None    
    inDataset1 = None        
    
if __name__ == '__main__':
    main()     