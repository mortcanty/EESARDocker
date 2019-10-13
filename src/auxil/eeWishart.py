'''
Created on 09.01.2017

The sequential omnibus algorithm for polarimetric SAR imagery

@author: mort
'''


import ee

def chi2cdf(chi2,df):
    ''' Chi square cumulative distribution function '''
    return ee.Image(chi2.divide(2)).gammainc(ee.Number(df).divide(2))

def det(image):
    '''return determinant of 1, 2, 4, or 9-band polarimetric image ''' 
    detmap = {'k':image.select(0),'ar':image.select(1),'ai':image.select(2),'pr':image.select(3),'pi':image.select(4),
              's':image.select(5),'br':image.select(6),'bi':image.select(7),'z':image.select(8)}
    expr = 'k*s*z+2*(ar*br*pr-ai*bi*pr+ai*br*pi+ar*bi*pi)-s*(pr*pr+pi*pi)-k*(br*br+bi*bi)-s*(ar*ar+ai*ai)'
    bands = image.bandNames().length()
    result = ee.Image(ee.Algorithms.If(bands.eq(1),image,None))
    result = ee.Image(ee.Algorithms.If(bands.eq(2),image.expression('b(0)*b(1)'),result))
    result = ee.Image(ee.Algorithms.If(bands.eq(4),image.expression('b(0)*b(3)-b(1)*b(1)-b(2)*b(2)'),result))
    return ee.Image(ee.Algorithms.If(bands.eq(9),image.expression(expr,detmap),result))

def log_det_sum(imList,j):
    '''return the log of the the determinant of the sum of the first j images in imList'''
    imList = ee.List(imList)
    sumj = ee.ImageCollection(imList.slice(0,j)).reduce(ee.Reducer.sum())                
    return ee.Image(det(sumj)).log()
    
def log_det(imList,j):
    '''return the log of the the determinant of the jth image in imList'''
    im = ee.Image(ee.List(imList).get(j.subtract(1)))
    return ee.Image(det(im)).log()
    
def pv(imList,median,j,enl):
    ''' calculate -2log(R_ell,j) and return P-value'''
    imList = ee.List(imList)
    p2 = ee.Image(imList.get(0)).bandNames().length()
    p = ee.Number(ee.Algorithms.If(p2.eq(2).Or(p2.eq(3)),p2,p2.sqrt()))
    j = ee.Number(j)
    f = p2
    one = ee.Number(1.0)
   
    rhoj = ee.Number(ee.Algorithms.If( p2.eq(2),
#      1 - (1.+1./(j*(j-1)))/(6.*p1*n) where p1=1                                      
        one.subtract(one.add(one.divide(j.multiply(j.subtract(one)))).divide(6).divide(enl)),
#      1 - (2*p2-1.)*(1.+1./(j*(j-1.))/6*p*n               
        one.subtract( p2.multiply(2).subtract(one) \
                      .multiply(one.add(one.divide(j.multiply(j.subtract(one))))) \
                      .divide(p.multiply(6).multiply(enl))
                    ) ))
   
    omega2j = ee.Number(ee.Algorithms.If( p2.eq(2),
#      -(f/4.)*(1.-1./rhoj)**2                                    
        one.subtract(one.divide(rhoj)).pow(2.0).multiply(f.divide(-4.0)),        
#      (f/4)*(1-1./rhoj)**2 + (1./(24.*n**2))*p2(p2-1.)*(1.+(2.*j-1.)/(j**2*(j-1.)**2))/rh0j**2       
        f.multiply(one.subtract(one.divide(rhoj)).pow(2)).divide(4). \
        add( one.divide(enl.pow(2).multiply(24)) \
             .multiply(p2.multiply(p2.subtract(one))) \
             .multiply(one.add(j.multiply(2)).subtract(one)) \
             .divide(j.pow(2).multiply(j.subtract(one)).pow(2)) \
             .divide(rhoj.pow(2))  ) ))
    
#  Zj = -2*lnRj
    Zj = ee.Image(log_det_sum(imList,j.subtract(1))) \
                 .multiply(j.subtract(1)) \
                 .add(log_det(imList,j))  \
                 .add(p.multiply(j).multiply(ee.Number(j).log())) \
                 .subtract(p.multiply(j.subtract(1)).multiply(j.subtract(1).log())) \
                 .subtract(ee.Image(log_det_sum(imList,j)).multiply(j)) \
                 .multiply(-2).multiply(enl)
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
    enl = ee.Number(prev.get('enl'))
    imList = prev.get('imList')
    pvs = ee.List(prev.get('pvs'))
    Z = ee.Image(prev.get('Z')) 
    pval,Zj = pv(imList,median,j,enl)  
# Z = sum_j Zj = -2lnQ_ell  
    Z = Z.add(Zj)
    return ee.Dictionary({'median':median,'imList':imList,'enl':enl,'pvs':pvs.add(pval),'Z':Z})   

def ells_iter(current,prev):
    ell = ee.Number(current)
    prev = ee.Dictionary(prev)
    pv_arr = ee.List(prev.get('pv_arr'))
    k = ee.Number(prev.get('k'))
    enl = ee.Number(prev.get('enl'))
    median = prev.get('median')
    imList = ee.List(prev.get('imList'))
#  number of bands (degrees of freedom)
    p2 = ee.Image(imList.get(0)).bandNames().length()
    imList_ell = imList.slice(ell.subtract(1))
    js = ee.List.sequence(2,k.subtract(ell).add(1))
    first = ee.Dictionary({'median':median,'imList':imList_ell,'enl':enl,'pvs':ee.List([]),'Z':ee.Image.constant(0.0)})
    result = ee.Dictionary(js.iterate(js_iter,first))
#  list of P-values for R_ell,j, j = ell+1 ... k    
    pvs = ee.List(result.get('pvs'))
#  omnibus test statistic -2lnQ_ell = sum_j(-2lnR_j)    
    Z =  ee.Image(result.get('Z'))
#  degrees of freedom 
    f = k.subtract(ell).multiply(p2)
    p = p2.sqrt()
    one = ee.Number(1)
#  rho, w2 values 
    rho = ee.Number(ee.Algorithms.If( p2.eq(2),
                                      one.subtract(
                                          k.divide(enl).subtract(one.divide(enl.multiply(k))) \
                                          .divide(k.subtract(one).multiply(6))),
                                      one.subtract( 
                                          p2.multiply(2).subtract(one) \
                                         .multiply(k.subtract(one.divide(k))) \
                                         .divide(k.subtract(one).multiply(p).multiply(6).multiply(enl)) )))
    w2 = ee.Number(ee.Algorithms.If( p2.eq(2),
                                     k.subtract(one).multiply(one.subtract(one.divide(rho)).pow(2)).divide(2),
                                     p2.multiply(p2.subtract(one)) \
                                     .multiply(k.subtract(one.divide(k.pow(2)))) \
                                     .divide(rho.pow(2).multiply(24).multiply(enl.pow(2))) \
                                     .subtract(p2.multiply(k.subtract(one)).multiply(one.subtract(one.divide(rho)).pow(2).divide(4))) ))
    Z = Z.multiply(rho)  
    PvQ = ee.Image.constant(1.0).subtract(chi2cdf(Z,f).multiply(ee.Number(1).subtract(w2))).subtract(chi2cdf(Z,f.add(4)).multiply(w2))
    PvQ = ee.Algorithms.If(median, PvQ.focal_median(),PvQ) 
#  put at end of current sequence     
    pvs = pvs.add(PvQ)          
    return ee.Dictionary({'k':k,'median':median,'enl':enl,'imList':imList,'pv_arr':pv_arr.add(pvs)})

def filter_j(current,prev):
    pv = ee.Image(current)
    prev = ee.Dictionary(prev)
    pvQ = ee.Image(prev.get('pvQ'))
    ell = ee.Number(prev.get('ell'))
    cmap = ee.Image(prev.get('cmap'))
    smap = ee.Image(prev.get('smap'))
    fmap = ee.Image(prev.get('fmap'))
    bmap = ee.Image(prev.get('bmap'))
    significance = ee.Image(prev.get('significance'))    
    j = ee.Number(prev.get('j'))
    cmapj = cmap.multiply(0).add(ell.add(j).subtract(1))
    cmap1 = cmap.multiply(0).add(1)
    tst = pv.lt(significance).And(pvQ.lt(significance)).And(cmap.eq(ell.subtract(1)))
    cmap = cmap.where(tst,cmapj)
    fmap = fmap.where(tst,fmap.add(1))
    smap = ee.Algorithms.If(ell.eq(1),smap.where(tst,cmapj),smap)
    idx = ell.add(j).subtract(2)
    tmp = bmap.select(idx)
    bname = bmap.bandNames().get(idx)
    tmp = tmp.where(tst,cmap1)
    tmp = tmp.rename([bname])    
    bmap = bmap.addBands(tmp,[bname],True)    
    return ee.Dictionary({'ell':ell,'j':j.add(1),'significance':significance,'pvQ':pvQ,'cmap':cmap,
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
    cmap = prev.get('cmap')
    smap = prev.get('smap')
    fmap = prev.get('fmap')
    bmap = prev.get('bmap')
    first = ee.Dictionary({'ell':ell,'j':1, 'significance':significance,'pvQ':pvQ,'cmap':cmap,
                                                                                  'smap':smap,
                                                                                  'fmap':fmap,
                                                                                  'bmap':bmap})     
    result = ee.Dictionary(ee.List(pvs).iterate(filter_j,first))   
    return ee.Dictionary({'ell':ell.add(1),'significance':significance,'cmap':result.get('cmap'),
                                                                       'smap':result.get('smap'),
                                                                       'fmap':result.get('fmap'),
                                                                       'bmap':result.get('bmap')})

def dmap_iter(current,prev):
    '''post-process for directional change maps'''
    prev = ee.Dictionary(prev)
    j = ee.Number(prev.get('j'))
    image = ee.Image(current) 
    p = image.bandNames().length()  
    avimg = ee.Image(prev.get('avimg'))
    diff = image.subtract(avimg)    
#  positive/negative definiteness from pivots      
    posd = ee.Image(ee.Algorithms.If(p.eq(1),det(diff).gt(0),False))
    posd = ee.Image(ee.Algorithms.If(p.eq(2).Or(p.eq(4)),det(diff.select(0)).gt(0).And(det(diff).gt(0)),posd))
    posd = ee.Image(ee.Algorithms.If(p.eq(9),det(diff.select(0)).gt(0).And(det(diff.select(0,1,2,5)).gt(0)).And(det(diff).gt(0)),posd))
    negd = ee.Image(ee.Algorithms.If(p.eq(1),det(diff).lt(0),False))
    negd = ee.Image(ee.Algorithms.If(p.eq(2).Or(p.eq(4)),det(diff.select(0)).lt(0).And(det(diff).gt(0)),negd))
    negd = ee.Image(ee.Algorithms.If(p.eq(9),det(diff.select(0)).lt(0).And(det(diff.select(0,1,2,5)).gt(0)).And(det(diff).lt(0)),negd))
    bmap = ee.Image(prev.get('bmap'))
    bmapj = bmap.select(j)
    dmap1 = bmapj.multiply(0).add(1)
    dmap2 = bmapj.multiply(0).add(2)
    dmap3 = bmapj.multiply(0).add(3)
    bmapj = bmapj.where(bmapj,dmap3)
    bmapj = bmapj.where(bmapj.And(posd),dmap1)
    bmapj = bmapj.where(bmapj.And(negd),dmap2)  
    bmap = bmap.addBands(bmapj,overwrite=True)
#  provisional means
    r = ee.Image(prev.get('r')).add(1)
    avimg = avimg.add(image.subtract(avimg).divide(r))
#  reset average image and r array if change occurred
    avimg = avimg.where(bmapj,image)
    r = r.where(bmapj,1)
    return ee.Dictionary({'avimg':avimg,'bmap':bmap,'j':j.add(1),'r':r})

def omnibus(imList,significance=0.0001,enl=4.4,median=False):
    '''return change maps for sequential omnibus change algorithm''' 
    imList = ee.List(imList)  
    k = imList.length()  
#  pre-calculate p-value array    
    ells = ee.List.sequence(1,k.subtract(1))
    first = ee.Dictionary({'k':k,'median':median,'enl':enl,'imList':imList,'pv_arr':ee.List([])}) 
    result = ee.Dictionary(ells.iterate(ells_iter,first))
    pv_arr = ee.List(result.get('pv_arr'))                                           
#  filter p-values to generate cmap, smap, fmap and bmap
    cmap = ee.Image(imList.get(0)).select(0).multiply(0.0)
    smap = ee.Image(imList.get(0)).select(0).multiply(0.0)
    fmap = ee.Image(imList.get(0)).select(0).multiply(0.0)   
    bmap = ee.Image.constant(ee.List.repeat(0,k.subtract(1)))    
    significance = ee.Image.constant(significance)
    first = ee.Dictionary({'ell':1,'significance':significance,'cmap':cmap,'smap':smap,'fmap':fmap,'bmap':bmap})
    result = ee.Dictionary(pv_arr.iterate(filter_ell,first))       
#  post-process bmap for change direction
    bmap = ee.Image(result.get('bmap')) 
    r = ee.Image(cmap.multiply(0).add(1))
    first = ee.Dictionary({'avimg':imList.get(0),'bmap':bmap,'j':ee.Number(0),'r':r})  
    dmap = ee.Dictionary(imList.slice(1).iterate(dmap_iter,first)).get('bmap')     
    return result.set('bmap',ee.Image(dmap))

if __name__ == '__main__':
    pass