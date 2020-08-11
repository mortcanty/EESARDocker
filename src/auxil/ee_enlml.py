#!/usr/bin/env python
#******************************************************************************
#  Name:     ee_enlml.py
#  Purpose: 
#    Estimation of ENL for polSAR covariance images, GEE version
#    using ML method with full covariance matrix (quad, dual or single)
#    If image is diagonal-only, just use first band (C11)
#    Anfinsen et al. (2009) IEEE TGARS 47(11), 3795-3809    
#    
#  Usage:
#    from auxil.ee_enlml import enl
# 
# Copyright (c) 2020 Mort Canty

import ee
import auxil.lookup as lookup
import numpy as np

def enl_iter(current,prev):
    idx = ee.Number(current)
    prev = ee.Dictionary(prev)
    enlhist = ee.List(prev.get('enlhist'))
    scale = ee.Number(prev.get('scale'))
    arrimg = ee.Image(prev.get('arrimg'))
    band = arrimg.select(idx)
    counter = band.multiply(0)
    count = counter.where(band.lt(0),1).rename('c').reduceRegion(ee.Reducer.sum(),scale=scale).get('c')
    return ee.Dictionary({'arrimg':arrimg,'scale':scale,'enlhist':enlhist.set(idx,count)})
       
def enl(image,scale=10):
#  construct the determinant image     
    detmap = {'k':image.select(0),'ar':image.select(1),'ai':image.select(2),'pr':image.select(3),'pi':image.select(4),
              's':image.select(5),'br':image.select(6),'bi':image.select(7),'z':image.select(8)}
    expr = 'k*s*z+2*(ar*br*pr-ai*bi*pr+ai*br*pi+ar*bi*pi)-s*(pr*pr+pi*pi)-k*(br*br+bi*bi)-s*(ar*ar+ai*ai)'
    bands = image.bandNames().length()
    result = ee.Image(ee.Algorithms.If(bands.lte(3),image.select(0),None))
    result = ee.Image(ee.Algorithms.If(bands.eq(4),image.expression('b(0)*b(3)-b(1)*b(1)-b(2)*b(2)'),result))
    detimg = ee.Image(ee.Algorithms.If(bands.eq(9),image.expression(expr,detmap),result))    
#  7x7 window average of log of determinant image        
    avlogdetimg = detimg.log().reduceNeighborhood(ee.Reducer.mean(),ee.Kernel.square(3.5),optimization='window') 
#  log of 7x7 wíndow average of the determinant image    
    logavdetimg = detimg.reduceNeighborhood(ee.Reducer.mean(),ee.Kernel.square(3.5),optimization='window').log()  
#  add the lookuptable as 500 bands    
    bands = bands.getInfo()
    if bands==9:
        d = 2
    elif bands==4:
        d = 1
    else:
        d = 0   
    lu = lookup.table()
    diff = avlogdetimg.subtract(logavdetimg)
    arrimg = diff.add(ee.Image.constant(list(lu[:500,d])))
#  shift and multiply to locate zero crossings
    arrimg1 = diff.add(ee.Image.constant(list(np.roll(lu[:500,d],1)))) 
    arrimg = arrimg.multiply(arrimg1)  
#  accumulate enl histogram
    first = ee.Dictionary({'arrimg':arrimg,'scale':scale,'enlhist':ee.List(list(np.zeros(500)))})   
    lst = ee.List(list(range(500)))   
    result = ee.List(ee.Dictionary(lst.iterate(enl_iter,first)).get('enlhist'))
    idx = 10 + d*10
    return result.set(idx,0).set(0,0)
   
if __name__ == '__main__':
    pass
    