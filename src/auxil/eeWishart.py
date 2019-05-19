'''
Created on 09.01.2017

The sequential omnibus algorithm for Sentinel-1 dual pol, diagonal only imagery

@author: mort
'''

ENL = 4.4

import ee
from eeMad import chi2cdf

def multbyenl(image):
    return ee.Image(image).multiply(ENL)

def log_det_sum(imList,j):
    '''return the log of the the determinant of the sum of the first j images in imList'''
    imList = ee.List(imList)
    nbands = ee.Image(imList.get(0)).bandNames().length() 
    sumj = ee.ImageCollection(imList.slice(0,j)).reduce(ee.Reducer.sum())
    return ee.Algorithms.If( nbands.eq(2),                         
        sumj.expression('b(0)*b(1)').log(),
        sumj.log() )                    
    
def log_det(imList,j):
    '''return the log of the the determinant of the jth image in imList'''
    im = ee.Image(ee.List(imList).get(j.subtract(1)))
    nbands = im.bandNames().length()
    return ee.Algorithms.If(nbands.eq(2),  
        im.expression('b(0)*b(1)').log(),
        im.log() )
    
def pv(imList,p,median,j):
    ''' calculate -2log(R_ell,j) and return P-value '''
    imList = ee.List(imList)
    p = ee.Number(p)
    j = ee.Number(j)
    f = p
    one = ee.Number(1.0)
# 1 - (1. + 1./(j*(j-1)))/(6.*p*n)    
    rhoj = one.subtract(one.add(one.divide(j.multiply(j.subtract(one)))).divide(6*ENL))
# -(f/4.)*(1.-1./rhoj)**2'    
    omega2j = one.subtract(one.divide(rhoj)).pow(2.0).multiply(f.divide(-4.0))
# Zj = -2*lnRj
    Zj = ee.Image(log_det_sum(imList,j.subtract(1))) \
                 .multiply(j.subtract(1)) \
                 .add(log_det(imList,j))  \
                 .add(p.multiply(j).multiply(ee.Number(j).log())) \
                 .subtract(p.multiply(j.subtract(1)).multiply(j.subtract(1).log())) \
                 .subtract(ee.Image(log_det_sum(imList,j)).multiply(j)) \
                 .multiply(-2*ENL)
# (1.-omega2j)*stats.chi2.cdf(rhoj*Zj,[f])+omega2j*stats.chi2.cdf(rhoj*Zj,[f+4])                 
    P = chi2cdf(Zj.multiply(rhoj),f).multiply(one.subtract(omega2j)) \
                 .add(chi2cdf(Zj.multiply(rhoj),f.add(4)).multiply(omega2j))
    PV = ee.Image.constant(1.0).subtract(P)
# 3x3 median filter    
    return (ee.Algorithms.If(median, PV.focal_median(), PV),Zj)    

def js_iter(current,prev):
    j = ee.Number(current)
    prev = ee.Dictionary(prev)
    median = prev.get('median')
    p = prev.get('p')
    imList = prev.get('imList')
    pvs = ee.List(prev.get('pvs'))
    Z = ee.Image(prev.get('Z')) 
    pval,Zj = pv(imList,p,median,j)  
# Z = sum_j Zj = -2lnQ_ell  
    Z = Z.add(Zj)
    return ee.Dictionary({'median':median,'p':p,'imList':imList,'pvs':pvs.add(pval),'Z':Z})   

def ells_iter(current,prev):
    ell = ee.Number(current)
    prev = ee.Dictionary(prev)
    pv_arr = ee.List(prev.get('pv_arr'))
    k = ee.Number(prev.get('k'))
    median = prev.get('median')
    p = prev.get('p')
    imList = ee.List(prev.get('imList'))
    imList_ell = imList.slice(ell.subtract(1))
    js = ee.List.sequence(2,k.subtract(ell).add(1))
    first = ee.Dictionary({'median':median,'p':p,'imList':imList_ell,'pvs':ee.List([]),'Z':ee.Image.constant(0.0)})
    result = ee.Dictionary(js.iterate(js_iter,first))
#  list of P-values for R_ell,j, j = ell+1 ... k    
    pvs = ee.List(result.get('pvs'))
#  omnibus test statistic -2lnQ_ell = sum_j(-2lnR_j)    
    Z =  ee.Image(result.get('Z'))
    f = k.subtract(ell).multiply(ee.Number(p))
    PvQ = ee.Image.constant(1.0).subtract(chi2cdf(Z,f)) 
    PvQ = ee.Algorithms.If(median, PvQ.focal_median(),PvQ)  
    pvs = pvs.add(PvQ)          
    return ee.Dictionary({'k':k,'p':p,'median':median,'imList':imList,'pv_arr':pv_arr.add(pvs)})

def filter_j(current,prev):
    pv = ee.Image(current)
    prev = ee.Dictionary(prev)
    pvQ = ee.Image(prev.get('pvQ'))
    ell = ee.Number(prev.get('ell'))
    cmap = ee.Image(prev.get('cmap'))
    smap = ee.Image(prev.get('smap'))
    fmap = ee.Image(prev.get('fmap'))
    bmap = ee.Image(prev.get('bmap'))
    dmaps = ee.List(prev.get('dmaps'))
    significance = ee.Image(prev.get('significance'))    
    j = ee.Number(prev.get('j'))
    
    dmap = ee.Image(dmaps.get(ell.add(j).subtract(2))) # ell+j-1 th element
    
    cmapj = cmap.multiply(0).add(ell.add(j).subtract(1))
#    cmap1 = cmap.multiply(0).add(1)
    tst = pv.lt(significance).And(pvQ.lt(significance)).And(cmap.eq(ell.subtract(1)))
    cmap = cmap.where(tst,cmapj)
    fmap = fmap.where(tst,fmap.add(1))
    smap = ee.Algorithms.If(ell.eq(1),smap.where(tst,cmapj),smap)
    idx = ell.add(j).subtract(2)
    tmp = bmap.select(idx)
    bname = bmap.bandNames().get(idx)
    tmp = tmp.where(tst,dmap)
    tmp = tmp.rename([bname])    
    bmap = bmap.addBands(tmp,[bname],True)    
    return ee.Dictionary({'ell':ell,'j':j.add(1),'significance':significance,'pvQ':pvQ,'dmaps':dmaps,
                                                                                       'cmap':cmap,
                                                                                       'smap':smap,
                                                                                       'fmap':fmap,
                                                                                       'bmap':bmap})

def filter_ell(current,prev):
    current = ee.List(current)
    pvs = current.slice(0,-1 )
    pvQ = ee.Image(current.get(-1))
    prev = ee.Dictionary(prev)
    ell = ee.Number(prev.get('ell'))
    significance = ee.Image(prev.get('significance'))
    useQ = ee.Number(prev.get('useQ'))
    pvQ = ee.Algorithms.If(useQ,pvQ,ee.Image.constant(0))
    dmaps = prev.get('dmaps')
    cmap = prev.get('cmap')
    smap = prev.get('smap')
    fmap = prev.get('fmap')
    bmap = prev.get('bmap')
    first = ee.Dictionary({'ell':ell,'j':1, 'significance':significance,'pvQ':pvQ,'dmaps':dmaps,
                                                                                  'cmap':cmap,
                                                                                  'smap':smap,
                                                                                  'fmap':fmap,
                                                                                  'bmap':bmap})     
    result = ee.Dictionary(ee.List(pvs).iterate(filter_j,first))   
    return ee.Dictionary({'ell':ell.add(1),'significance':significance,'useQ':useQ,'dmaps':dmaps,
                                                                                   'cmap':result.get('cmap'),
                                                                                   'smap':result.get('smap'),
                                                                                   'fmap':result.get('fmap'),
                                                                                   'bmap':result.get('bmap')})

def dmap_iter(current,prev):
    im1 = ee.Image(current)
    prev = ee.Dictionary(prev)
    dmaps = ee.List(prev.get('dmaps'))
    im2list = ee.List(prev.get('im2list'))
    im2 = ee.Image(im2list.get(0))
    eiv1 = im2.subtract(im1).select(0)
    eiv2 = im2.subtract(im1).select(1)
    dmap = ee.Image(im1.select(0)).multiply(0.0).rename(['direction'])
    dmap1 = dmap.add(1)
    dmap2 = dmap.add(2)
    dmap = dmap.add(3)
    tst1 = eiv1.gt(0).And(eiv2.gt(0)) 
    tst2 = eiv1.lt(0).And(eiv2.lt(0)) 
    dmap = dmap.where(tst1,dmap1)
    dmap = dmap.where(tst2,dmap2)
    return ee.Dictionary({'im2list':im2list.remove(im2),'dmaps':dmaps.add(dmap)})

def omnibus(imList,significance=0.0001,median=False,useQ=False):
    '''return change maps for sequential omnibus change algorithm''' 
    imList = ee.List(imList)  
    p = ee.Image(imList.get(0)).bandNames().length()
    k = imList.length()  
#  generate change direction map    
    first = ee.Dictionary({'im2list':imList.slice(1),'dmaps':ee.List([])})
    result = ee.Dictionary(imList.slice(0,-1).iterate(dmap_iter,first))
    dmaps = ee.List(result.get('dmaps'))   
#  pre-calculate p-value array    
    imList = imList.map(multbyenl)   
    ells = ee.List.sequence(1,k.subtract(1))
    first = ee.Dictionary({'k':k,'p':p,'median':median,'imList':imList,'pv_arr':ee.List([])})
    result = ee.Dictionary(ells.iterate(ells_iter,first))
    pv_arr = ee.List(result.get('pv_arr'))  
#  filter p-values to generate cmap, smap, fmap and bmap
    cmap = ee.Image(imList.get(0)).select(0).multiply(0.0)
    smap = ee.Image(imList.get(0)).select(0).multiply(0.0)
    fmap = ee.Image(imList.get(0)).select(0).multiply(0.0)   
    bmap = ee.Image.constant(ee.List.repeat(0,k.subtract(1)))    
    significance = ee.Image.constant(significance)
    first = ee.Dictionary({'ell':1,'significance':significance,'useQ':useQ,'dmaps':dmaps,
                                                                           'cmap':cmap,
                                                                           'smap':smap,
                                                                           'fmap':fmap,
                                                                           'bmap':bmap})
    return ee.Dictionary(pv_arr.iterate(filter_ell,first)) 

if __name__ == '__main__':
    pass