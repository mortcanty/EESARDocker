#!/usr/bin/env python
#******************************************************************************
#  Name:     plr.py
#  Purpose:  Module for probabilistic label relaxation 
#  Usage:             
#    python plr.py [options] class_prob_image 
#
#  Copyright (c) 2018, Mort Canty

import os, time, sys, getopt
import numpy as np
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly, GDT_Byte
    
def plr(infile,nitr=3):    
    path = os.path.dirname(infile)
    basename = os.path.basename(infile)
    root, ext = os.path.splitext(basename)
    outfile = path+'/'+root+'_plr'+ext  
    gdal.AllRegister()                  
    inDataset = gdal.Open(infile,GA_ReadOnly)     
    cols = inDataset.RasterXSize
    rows = inDataset.RasterYSize    
    classes = inDataset.RasterCount
    print('=====================')
    print('       PLR')
    print('=====================')
    print('infile:  %s'%infile)
    print('iterations:  %i'%nitr)
    start = time.time()                                   
    prob_image = np.zeros((classes,rows,cols))
    for k in range (classes):
        band = inDataset.GetRasterBand(k+1)
        prob_image[k,:,:] = np.array(band.ReadAsArray(0,0,cols,rows),dtype=np.float64)/255.
#  compatibility matrix
    Pmn = np.zeros((classes,classes))
    n_samples = (cols-1)*(rows-1)
    samplem = np.reshape(prob_image[:,0:rows-1,0:cols-1],(classes,n_samples))
    samplen = np.reshape(prob_image[:,1:rows,0:cols-1],(classes,n_samples))
    sampleu = np.reshape(prob_image[:,0:rows-1,1:cols],(classes,n_samples))
    print('estimating compatibility matrix...')
    for j in range(n_samples):     
        m1 = np.argmax(samplem[:,j])
        n1 = np.argmax(samplen[:,j])
        u1 = np.argmax(sampleu[:,j])
        Pmn[m1,n1] += 1  
        Pmn[m1,u1] += 1  
    for j in range(classes):
        n = np.sum(Pmn[j,:])
        if n>0:
            Pmn[j,:] /= n   
    itr = 0
    temp = prob_image*0
    print('label relaxation...')
    while itr<nitr:
        print('iteration %i'%(itr+1))
        Pm = np.zeros(classes)
        Pn = np.zeros(classes)
        for i in range(1,rows-1):
            for j in range(1,cols-1):
                Pm[:] = prob_image[:,i,j]
                Pn[:] = prob_image[:,i-1,j]/4
                Pn[:] += prob_image[:,i+1,j]/4
                Pn[:] += prob_image[:,i,j-1]/4
                Pn[:] += prob_image[:,i,j+1]/4
                Pn = np.transpose(Pn)
                den = (np.dot(np.dot(Pm,Pmn),Pn))
                if den == 0:
                    Pm_new = Pm
                else:
                    Pm_new = Pm*(np.dot(Pmn,Pn))/den
                temp[:,i,j] = Pm_new
        prob_image = temp
        itr += 1     
    class_image = np.argmax(prob_image,axis=0)+1   
#  write to disk
    driver = gdal.GetDriverByName('GTiff')    
    outDataset = driver.Create(outfile,cols,rows,1,GDT_Byte)
    projection = inDataset.GetProjection()
    geotransform = inDataset.GetGeoTransform()
    if geotransform is not None:
        outDataset.SetGeoTransform(geotransform)
    if projection is not None:
        outDataset.SetProjection(projection)               
    outBand = outDataset.GetRasterBand(1)
    outBand.WriteArray(class_image,0,0) 
    outBand.FlushCache() 
    outDataset = None
    inDataset = None
    print('result written to: '+outfile)    
    print('elapsed time: '+str(time.time()-start))                        
    print('--done------------------------')      
    
def main():
    usage = '''Usage: python %s [-h] [-i iterations] probfileName
            '''%sys.argv[0]
    usage = '''
Usage:
------------------------------------------------

Probabilistic label relaxation post processing  

python %s [OPTIONS]  classProbFileName

Options:
  
  -h         this help  
  -i  <int>  number of iterations (default 3)

-------------------------------------------------'''%sys.argv[0]                  
    options,args = getopt.getopt(sys.argv[1:],'hi:')
    iterations = 3
    for option, value in options: 
        if option == '-h':
            print(usage)
            return 
        elif option == '-i':
            iterations = eval(value)  
    infile = args[0] 
    plr(infile,iterations)
              
if __name__ == '__main__':
    main()    