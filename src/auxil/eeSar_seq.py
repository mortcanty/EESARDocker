'''
Created on 21.06.2018

@author: mort

ipywidget interface to the GEE for sequential SAR change detection

'''
import ee, time, warnings, math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import ipywidgets as widgets
from IPython.display import display
from ipyleaflet import (Map,DrawControl,TileLayer,
                        basemaps,basemap_to_tiles,
                        LayersControl,
                        MeasureControl,
                        FullScreenControl)
from auxil.eeWishart import omnibus
from auxil.eeRL import refinedLee
from auxil.ee_enlml import enl
from geopy.geocoders import photon

ee.Initialize()

warnings.filterwarnings("ignore")

poly = ee.Geometry.MultiPolygon([])

def update_figure(fig,line,profile):    
    fig.title = 'Change Profile'
    line.x = range(1,count)
    line.y = profile

geolocator = photon.Photon(timeout=10)

def get_incidence_angle(image):
    ''' grab the mean incidence angle '''
    result = ee.Image(image).select('angle') \
           .reduceRegion(ee.Reducer.mean(),geometry=poly,maxPixels=1e9) \
           .get('angle') \
           .getInfo()
    if result is not None:
        return round(result,2)
    else:
#      incomplete overlap, so use all of the image geometry        
        return round(ee.Image(image).select('angle') \
           .reduceRegion(ee.Reducer.mean(),maxPixels=1e9) \
           .get('angle') \
           .getInfo(),2)

def get_vvvh(image):   
    ''' get 'VV' and 'VH' bands from sentinel-1 imageCollection and restore linear signal from db-values '''
    return image.select('VV','VH').multiply(ee.Image.constant(math.log(10.0)/10.0)).exp()

def get_image(current,image):
    ''' accumulate a single image from a collection of images '''
    return ee.Image.cat(ee.Image(image),current)    
    
def clipList(current,prev):
    ''' clip a list of images and multiply by ENL'''
    imlist = ee.List(ee.Dictionary(prev).get('imlist'))
    poly = ee.Dictionary(prev).get('poly') 
    enl = ee.Number(ee.Dictionary(prev).get('enl')) 
    ctr = ee.Number(ee.Dictionary(prev).get('ctr'))   
    stride = ee.Number(ee.Dictionary(prev).get('stride'))
    imlist =  ee.Algorithms.If(ctr.mod(stride).eq(0),
        imlist.add(ee.Image(current).multiply(enl).clip(poly)),imlist)
    return ee.Dictionary({'imlist':imlist,'poly':poly,'enl':enl,'ctr':ctr.add(1),'stride':stride})

def makefeature(data):
    ''' for exporting as CSV to Drive '''
    return ee.Feature(None, {'data': data})

def handle_draw(self, action, geo_json):
    global poly
    coords =  geo_json['geometry']['coordinates']
    if action == 'created':
        poly = ee.Geometry.MultiPolygon(poly.coordinates().add(coords))
        w_preview.disabled = True
        w_export_ass.disabled = True
        w_collect.disabled = False
        w_ENL.disabled = False
    elif action == 'deleted':
        poly1 = ee.Geometry.MultiPolygon(coords)
        poly = poly.difference(poly1)
        if len(poly.coordinates().getInfo()) == 0:
            w_collect.disabled = True            

def GetTileLayerUrl(ee_image_object):
    map_id = ee.Image(ee_image_object).getMapId()
    return map_id["tile_fetcher"].url_format

w_collection = widgets.Text(
    value='COPERNICUS/S1_GRD',
    placeholder=' ',
    description='Collection:',
    disabled=False
)
w_enl = widgets.BoundedFloatText(
    value='4.4',
    min=3.0,
    max=20.0,
    step=0.1,
    description='ENL:',
    disabled=False
)
w_location = widgets.Text(
    value='Jülich',
    placeholder=' ',
    description='',
    disabled=False
)
w_orbitpass = widgets.RadioButtons(
    options=['ASCENDING','DESCENDING'],
    value='ASCENDING',
    description='Orbit pass:',
    disabled=False
)
w_changemap = widgets.RadioButtons(
    options=['Bitemporal','First','Last','Frequency'],
     value='First',
    description='Map:',
    disabled=False
)
w_bmap = widgets.BoundedIntText(
    layout = widgets.Layout(width='50px'),
    min=1,
    value=1,
    description='',
    disabled=True
)
w_platform = widgets.RadioButtons(
    layout = widgets.Layout(width='25%'),
    options=['Both','A','B'],
     value='Both',
    description='Platform:',
    disabled=False
)
w_relativeorbitnumber = widgets.IntText(
    layout = widgets.Layout(width='200px'),
    value='0',
    description='Rel orbit:',
    disabled=False
)
w_exportassetsname = widgets.Text(
    layout = widgets.Layout(width='200px'),
    value='users/<username>/<path>',
    placeholder=' ',
    disabled=False
)
w_exportdrivename = widgets.Text(
    layout = widgets.Layout(width='200px'),
    value='gee/<path>',
    placeholder=' ',
    disabled=False
)
w_exportscale = widgets.FloatText(
    layout = widgets.Layout(width='150px'),
    value=10,
    placeholder=' ',
    description='Scale ',
    disabled=False
)
w_startdate = widgets.Text(
    value='2018-04-01',
    placeholder=' ',
    description='Start date:',
    disabled=False
)
w_enddate = widgets.Text(
    value='2018-11-01',
    placeholder=' ',
    description='End date:',
    disabled=False
)
w_stride = widgets.BoundedIntText(
    layout = widgets.Layout(width='200px'),
    value=1,
    min=1,
    description='Stride:',
    disabled=False
)
w_median = widgets.Checkbox(
    value=True,
    description='5x5 Median filter',
    disabled=False
)
w_S2 = widgets.Checkbox(
    layout = widgets.Layout(width='200px'),
    value=False,
    description='Show best S2',
    disabled=False
)
w_Q = widgets.Checkbox(
    layout = widgets.Layout(width='200px'),
    value=True,
    description='Quick Preview',
    disabled=False
)
w_significance = widgets.BoundedFloatText(
    layout = widgets.Layout(width='200px'),
    value='0.01',
    min=0.001,
    max=0.05,
    step=0.001,
    description='Significance:',
    disabled=False
)
w_opacity = widgets.BoundedFloatText(
    layout = widgets.Layout(width='200px'),
    value='1.0',
    min=0.0,
    max=1.0,
    step=0.1,
    description='Opacity:',
    disabled=False
)
w_maskchange = widgets.Checkbox(
    value=False,
    description='NC mask',
    disabled=False
)
w_maskwater = widgets.Checkbox(
    value=True,
    description='Water mask',
    disabled=False
)
w_out = widgets.Output(
    layout={'border': '1px solid black'}
)

w_drive = widgets.Text(
    value='<path>',
    placeholder=' ',
    description='',
    disabled=False
)

w_goto = widgets.Button(description='GoTo',disabled=False)
w_export_atsf = widgets.Button(description='ATSFtoDrive',disabled=True)
w_export_s2 = widgets.Button(description='S2toDrive',disabled=True)
w_coll = widgets.HBox([w_collection,w_enl,w_goto,w_location])
w_opac = widgets.VBox([w_opacity,w_maskchange,w_maskwater],layout = widgets.Layout(width='50%'))
w_collect = widgets.Button(description="Collect",disabled=True)
w_preview = widgets.Button(description="Preview",disabled=True)
w_review = widgets.Button(description="Review",disabled=False)
w_plot = widgets.Button(description='PlotFromAsset',disabled=False)
w_clearpoly = widgets.Button(description="ClearPoly")
w_ENL = widgets.Button(description="EstimateENL",disabled=True)
w_export_ass = widgets.Button(description='ChangesToAssets',disabled=True)
w_export_drv = widgets.Button(description='ChangesToDrive',disabled=True)
w_dates = widgets.HBox([w_relativeorbitnumber,w_startdate,w_enddate,w_stride])
w_change = widgets.HBox([w_changemap,w_bmap])
w_orbit = widgets.HBox([w_orbitpass,w_platform,w_change,w_opac])
w_signif = widgets.HBox([w_significance,w_S2,w_Q,w_median,w_exportscale],layout = widgets.Layout(width='99%'))
w_run = widgets.HBox([w_collect,w_preview,w_plot,w_clearpoly,w_ENL,w_review])
w_reset = widgets.Button(description='Reset',disabled=False)
w_output = widgets.HBox([w_reset,w_out])
w_exp = widgets.HBox([w_export_ass,w_exportassetsname,w_export_drv,w_drive,w_export_atsf,w_export_s2])


box = widgets.VBox([w_output,w_coll,w_dates,w_orbit,w_signif,w_run,w_exp])

def on_widget_change(b):
    w_preview.disabled = True
    w_export_ass.disabled = True
    w_export_drv.disabled = True
    w_export_atsf.disabled = True
    
def on_changemap_widget_change(b):   
    if b['new']=='Bitemporal':
        w_bmap.disabled=False
    else:
        w_bmap.disabled=True
    
w_orbitpass.observe(on_widget_change,names='value')
w_platform.observe(on_widget_change,names='value')
w_relativeorbitnumber.observe(on_widget_change,names='value')
w_startdate.observe(on_widget_change,names='value')
w_enddate.observe(on_widget_change,names='value')
w_stride.observe(on_widget_change,names='value')
w_collection.observe(on_widget_change,names='value')
w_enl.observe(on_widget_change,names='value')
w_enddate.observe(on_widget_change,names='value')
w_median.observe(on_widget_change,names='value')
w_significance.observe(on_widget_change,names='value')
w_changemap.observe(on_changemap_widget_change,names='value')


def getS1collection():
    s1 =  ee.ImageCollection('COPERNICUS/S1_GRD') \
                      .filterBounds(poly) \
                      .filterDate(ee.Date(w_startdate.value), ee.Date(w_enddate.value)) \
                      .filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV','VH'])) \
                      .filter(ee.Filter.eq('resolution_meters', 10)) \
                      .filter(ee.Filter.eq('instrumentMode', 'IW')) \
                      .filter(ee.Filter.eq('orbitProperties_pass', w_orbitpass.value))    
    return s1.filter(ee.Filter.contains(rightValue=poly,leftField='.geo'))

def getS2collection():
    s2 = ee.ImageCollection('COPERNICUS/S2') \
                      .filterBounds(poly) \
                      .filterDate(ee.Date(w_startdate.value),ee.Date(w_enddate.value)) \
                      .sort('CLOUDY_PIXEL_PERCENTAGE',True)
    return s2.filter(ee.Filter.contains(rightValue=poly,leftField='.geo'))                      
                      
def on_reset_button_clicked(b):
    with w_out:
        w_out.clear_output()
        print('Algorithm output')   
        
w_reset.on_click(on_reset_button_clicked)                           

# def on_poly_button_clicked(b):
#     asset = w_exportassetsname.value
#     assetpoly = ee.Image(asset).select(0).geometry()
#     with w_out:
#         w_out.clear_output()
#         print(str(assetpoly.bounds().getInfo()))
        

def on_clearpoly_button_clicked(b):
    global poly
    poly = ee.Geometry.MultiPolygon([])
    with w_out:
        w_out.clear_output()
        print('Algorithm output')    
    w_collect.disabled = True
    w_preview.disabled = True
    w_ENL.disabled = True      
    
w_clearpoly.on_click(on_clearpoly_button_clicked)    
 

def on_ENL_button_clicked(b):
    with w_out:
        try:
            w_out.clear_output()            
            print('ENL calculation for %s ...'%timestamplist1[0])        
            y = np.array(enl(collectionfirst.clip(poly),w_exportscale.value).getInfo())                     
            x = np.linspace(0,50,500)
            itp = interp1d(x,y, kind='linear')
            window_size, poly_order = 21, 3
            y_sg = savgol_filter(itp(x), window_size, poly_order) 
            idx = np.where(y_sg<0)  
            y_sg[idx] = 0         
            print('mode: %s'%str(np.argmax(y_sg)/10.0))        
            ax = plt.subplot(111)
            ax.plot(x,y,label = 'ENL')
            ax.plot(x,y_sg,label = 'ENL smoothed')
            ax.legend()   
            ax.grid() 
            plt.show()           
        except Exception as e:
            print('Error: %s'%e)     
    
w_ENL.on_click(on_ENL_button_clicked)    

def on_collect_button_clicked(b):
    global result,collection,count,imList,poly,timestamplist1,timestamps2, \
           s2_image,rons,mean_incidence,collectionmosaic,collectionfirst,archive_crs,coords,wc 
    with w_out:
        try:
            if (w_collection.value == 'COPERNICUS/S1_GRD') or (w_collection.value == ''): 
                w_out.clear_output()
                print('running on GEE archive COPERNICUS/S1_GRD (please wait for raster overlay) ...')
#               coords = ee.List(poly.bounds().coordinates().get(0))
                collection = getS1collection()              
                if w_relativeorbitnumber.value > 0:
                    collection = collection.filter(ee.Filter.eq('relativeOrbitNumber_start', int(w_relativeorbitnumber.value)))   
                if w_platform.value != 'Both':
                    collection = collection.filter(ee.Filter.eq('platform_number', w_platform.value))         
                collection = collection.sort('system:time_start') 
                acquisition_times = ee.List(collection.aggregate_array('system:time_start')).getInfo()              
                count = len(acquisition_times)      
                if count<2:
                    raise ValueError('Less than 2 images found')
                timestamplist = []
                for timestamp in acquisition_times:
                    tmp = time.gmtime(int(timestamp)/1000)
                    timestamplist.append(time.strftime('%x', tmp))  
#              make timestamps in YYYYMMDD format            
                timestamplist = [x.replace('/','') for x in timestamplist]  
                timestamplist = ['T20'+x[4:]+x[0:4] for x in timestamplist]         
                timestamplist = timestamplist[::int(w_stride.value)]
#              in case of duplicates add running integer
                timestamplist1 = [timestamplist[i] + '_' + str(i+1) for i in range(len(timestamplist))]     
                count = len(timestamplist)
                if count<2:
                    raise ValueError('Less than 2 images found, decrease stride')            
                relativeorbitnumbers = map(int,ee.List(collection.aggregate_array('relativeOrbitNumber_start')).getInfo())
                rons = list(set(relativeorbitnumbers))
                print('Images found: %i, platform: %s'%(count,w_platform.value))
                print('Number of 10m pixels contained: %i'%math.floor(poly.area().getInfo()/100.0))
                print('Acquisition dates: %s to %s'%(str(timestamplist[0]),str(timestamplist[-1])))
                print('Relative orbit numbers: '+str(rons))
                if len(rons)==1:
                    mean_incidence = get_incidence_angle(collection.first())
                    print('Mean incidence angle: %f'%mean_incidence)
                else:
                    mean_incidence = 'undefined'
                    print('Mean incidence angle: (select one rel. orbit)')
                pcollection = collection.map(get_vvvh)            
                collectionfirst = ee.Image(pcollection.first())
                w_exportscale.value = collectionfirst.projection().nominalScale().getInfo()          
                pList = pcollection.toList(500)   
                first = ee.Dictionary({'imlist':ee.List([]),'poly':poly,'enl':ee.Number(w_enl.value),'ctr':ee.Number(0),'stride':ee.Number(int(w_stride.value))}) 
                imList = ee.List(ee.Dictionary(pList.iterate(clipList,first)).get('imlist'))              
#              get a vorschau as collection mean                                           
                collectionmosaic = collection.mosaic().select(0,1).rename('b0','b1')
                percentiles = collectionmosaic.reduceRegion(ee.Reducer.percentile([2,98]),geometry=poly,scale=w_exportscale.value,maxPixels=10e9)
                mn = ee.Number(percentiles.get('b0_p2'))
                mx = ee.Number(percentiles.get('b0_p98'))        
                vorschau = collectionmosaic.select(0).visualize(min=mn, max=mx, opacity=w_opacity.value) 
            else:
                w_out.clear_output()
                collection = ee.ImageCollection(w_collection.value)
                print('running on local collection %s \n ignoring start and end dates (please wait for raster overlay) ...'%w_collection.value)  
                count = collection.size().getInfo()  
                print('Images found: %i'%count )          
                collectionfirst = ee.Image(collection.first())  
                poly = collectionfirst.geometry()   
#                coords = ee.List(poly.bounds().coordinates().get(0))   
                center = poly.centroid().coordinates().getInfo()
                center.reverse()
                m.center = center                
                w_exportscale.value = collectionfirst.projection().nominalScale().getInfo()
                if collectionfirst.get('system:time_start').getInfo() is not None:
                    acquisition_times = ee.List(collection.aggregate_array('system:time_start')).getInfo()  
                    timestamplist1 = []
                    for timestamp in acquisition_times:
                        tmp = time.gmtime(int(timestamp)/1000)
                        timestamplist1.append(time.strftime('%x', tmp))            
                    timestamplist1 = [x.replace('/','') for x in timestamplist1]  
                    timestamplist1 = ['T20'+x[4:]+x[0:4] for x in timestamplist1]    
                    print('Acquisition dates: %s'%str(timestamplist1))    
                else:
                    timestamplist1 = ['T%i'%(i+1) for i in range(count)]
                    print('No time property available: acquisitions: %s'%str(timestamplist1))         
#              get a vorschau from collection mean                 
                collectionmosaic = collection.mosaic().clip(poly)
                percentiles = collectionmosaic.select(0).rename('b0').reduceRegion(ee.Reducer.percentile([2,98]),scale=w_exportscale.value,maxPixels=10e9)
                mn = ee.Number(percentiles.get('b0_p2'))
                mx = ee.Number(percentiles.get('b0_p98'))        
                vorschau = collectionmosaic.select(0).visualize(min=mn, max=mx, opacity=w_opacity.value)       
                imList = collection.toList(100)
#          get GEE S1 archive crs for eventual image series export               
#            archive_crs = ee.Image(getS1collection(coords).first()).select(0).projection().crs().getInfo()
            archive_crs = ee.Image(getS1collection().first()).select(0).projection().crs().getInfo()
#          run the algorithm        
            result = omnibus(imList,w_significance.value,w_enl.value,w_median.value)         
            w_preview.disabled = False
            w_ENL.disabled = False
            w_export_atsf.disabled = True
            s2_image = None
#          display collection or S2 
            if len(m.layers)>3:
                m.remove_layer(m.layers[3])
            if w_S2.value:
#              display sentinel-2 if available              
                collection2 = getS2collection() 
                count1 = collection2.size().getInfo()
                if count1>0:    
                    s2_image =  ee.Image(collection2.first()).select(['B2','B3','B4'])      
                    percentiles = s2_image.reduceRegion(ee.Reducer.percentile([2,98]),scale=w_exportscale.value,maxPixels=10e9)         
                    mn = percentiles.values(['B2_p2','B3_p2','B4_p2'])
                    mx = percentiles.values(['B2_p98','B3_p98','B4_p98'])
                    vorschau = s2_image.visualize(min=mn,max=mx,opacity=w_opacity.value)           
                    timestamp = s2_image.get('system:time_start').getInfo() 
                    timestamp = time.gmtime(int(timestamp)/1000)
                    timestamp = time.strftime('%x', timestamp).replace('/','')
                    timestamps2 = '20'+timestamp[4:]+timestamp[0:4]
                    print('Sentinel-2 from %s'%timestamps2) 
                    w_export_s2.disabled = False
            m.add_layer(TileLayer(url=GetTileLayerUrl(vorschau)))
          
        except Exception as e:
            print('Error: %s'%e)       

w_collect.on_click(on_collect_button_clicked)

def on_goto_button_clicked(b):
    try:
        location = geolocator.geocode(w_location.value)
        m.center = (location.latitude,location.longitude)
        m.zoom = 11
    except Exception as e:
        with w_out:
            print('Error: %s'%e)

w_goto.on_click(on_goto_button_clicked)

def on_preview_button_clicked(b):
    global cmap,smap,fmap,bmap,avimgs,avimg,pvQ,avimglog,count,watermask
    watermask = ee.Image('UMD/hansen/global_forest_change_2015').select('datamask').eq(1)  
    with w_out:  
        try:       
            jet = 'black,blue,cyan,yellow,red'
            rcy = 'black,red,cyan,yellow'
            smap = ee.Image(result.get('smap')).byte()
            cmap = ee.Image(result.get('cmap')).byte()
            fmap = ee.Image(result.get('fmap')).byte() 
            bmap = ee.Image(result.get('bmap')).byte()   
#          the atsf                    
            avimgs = ee.List(result.get('avimgs'))
            avimg = ee.Image(ee.List(result.get('avimgs')).get(-1)).clip(poly)  
            avimglog = ee.Image(result.get('avimglog')).byte().clip(poly)     
#          for control           
            pvQ =  ee.Image(result.get('pvQ'))              
            palette = jet
            w_out.clear_output()
            print('Series length: %i images, previewing (please wait for raster overlay) ...'%count)
            if w_changemap.value=='First':
                mp = smap
                mx = count
                print('Interval of first change:\n blue = early, red = late')
            elif w_changemap.value=='Last':
                mp=cmap
                mx = count
                print('Interval of last change:\n blue = early, red = late')
            elif w_changemap.value=='Frequency':
                mp = fmap
                mx = count/2
                print('Change frequency :\n blue = few, red = many')
            else:
                sel = int(w_bmap.value)
                sel = min(sel,count-1)
                sel = max(sel,1)
                print('Bitemporal: %s-->%s'%(timestamplist1[sel-1],timestamplist1[sel]))
                print('red = positive definite, cyan = negative definite, yellow = indefinite')     
                mp = bmap.select(sel-1).clip(poly)
                palette = rcy
                mx = 3     
            if len(m.layers)>3:
                m.remove_layer(m.layers[3])
            if not w_Q.value:
                mp = mp.reproject(crs=archive_crs,scale=float(w_exportscale.value))
            if w_maskwater.value==True:
                mp = mp.updateMask(watermask)
            if w_maskchange.value==True:    
                mp = mp.updateMask(mp.gt(0))    
            m.add_layer(TileLayer(url=GetTileLayerUrl(mp.visualize(min=0, max=mx, palette=palette,opacity = w_opacity.value))))
            w_export_ass.disabled = False
            w_export_drv.disabled = False
            w_export_atsf.disabled = False
        except Exception as e:
            print('Error: %s'%e)
    
w_preview.on_click(on_preview_button_clicked)   

def on_review_button_clicked(b):
    global poly
    watermask = ee.Image('UMD/hansen/global_forest_change_2015').select('datamask').eq(1)  
    with w_out:  
        try:       
            asset = ee.Image(w_exportassetsname.value)
            poly = ee.Geometry.Polygon(ee.Geometry(asset.get('system:footprint')).coordinates())
            center = poly.centroid().coordinates().getInfo()
            center.reverse()
            m.center = center  
            bnames = asset.bandNames().getInfo()[3:-2]
            count = len(bnames)               
            jet = 'black,blue,cyan,yellow,red'
            rcy = 'black,red,cyan,yellow'
            smap = asset.select('smap').byte()
            cmap = asset.select('cmap').byte()
            fmap = asset.select('fmap').byte()
            bmap = asset.select(list(range(3,count+3)),bnames).byte()      
            palette = jet
            w_out.clear_output()
            print('Series length: %i images, reviewing (please wait for raster overlay) ...'%(count+1))
            if w_changemap.value=='First':
                mp = smap
                mx = count
                print('Interval of first change:\n blue = early, red = late')
            elif w_changemap.value=='Last':
                mp = cmap
                mx = count
                print('Interval of last change:\n blue = early, red = late')
            elif w_changemap.value=='Frequency':
                mp = fmap
                mx = count/2
                print('Change frequency :\n blue = few, red = many')
            else:
                sel = int(w_bmap.value)-1
                sel = min(sel,count-1)
                sel = max(sel,0)
                if sel>0:
                    print('Bitemporal: %s --> %s'%(bnames[sel-1],bnames[sel]))
                else:
                    print('Bitemporal: image1 --> %s'%bnames[sel])
                print('red = positive definite, cyan = negative definite, yellow = indefinite')     
                mp = bmap.select(sel)
                palette = rcy
                mx = 3     
            if len(m.layers)>3:
                m.remove_layer(m.layers[3])
            if w_maskwater.value==True:
                mp = mp.updateMask(watermask)
            if w_maskchange.value==True:    
                mp = mp.updateMask(mp.gt(0))    
            m.add_layer(TileLayer(url=GetTileLayerUrl(mp.visualize(min=0, max=mx, palette=palette,opacity = w_opacity.value))))
            w_export_ass.disabled = False
            w_export_drv.disabled = False
            w_export_atsf.disabled = False
            w_collect.disabled = False
        except Exception as e:
            print('Error: %s'%e)
    
w_review.on_click(on_review_button_clicked)   

def on_plot_button_clicked(b):          
#  plot change fractions        
    global bmap1 
    watermask = ee.Image('UMD/hansen/global_forest_change_2015').select('datamask').eq(1) 
    def plot_iter(current,prev):
        current = ee.Image.constant(current)
        plots = ee.List(prev) 
        res = bmap1.multiply(0) \
                  .where(bmap1.eq(current),1) \
                  .reduceRegion(ee.Reducer.mean(),scale=w_exportscale.value,maxPixels=10e10)
        return ee.List(plots.add(res))
    with w_out:
        try:
            w_out.clear_output()            
            print('Change fraction plots ...')                  
            assetImage = ee.Image(w_exportassetsname.value)
            k = assetImage.bandNames().length().subtract(4).getInfo()            
            bmap1 = assetImage.select(ee.List.sequence(3,k+2)).updateMask(watermask)             
            plots = ee.List(ee.List([1,2,3]).iterate(plot_iter,ee.List([]))).getInfo()           
            bns = np.array(list([s[3:9] for s in list(plots[0].keys())])) 
            x = range(1,k+1)  
            _ = plt.figure(figsize=(10,5))
            plt.plot(x,list(plots[0].values()),'ro-',label='posdef')
            plt.plot(x,list(plots[1].values()),'co-',label='negdef')
            plt.plot(x,list(plots[2].values()),'yo-',label='indef')        
            ticks = range(0,k+2)
            labels = [str(i) for i in range(0,k+2)]
            labels[0] = ' '
            labels[-1] = ' '
            labels[1:-1] = bns 
            if k>50:
                for i in range(1,k+1,2):
                    labels[i] = ''
            plt.xticks(ticks,labels,rotation=90)
            plt.legend()
            fn = w_exportassetsname.value.replace('/','-')+'.png'
            plt.savefig(fn,bbox_inches='tight') 
            w_out.clear_output()
            plt.show()
            print('Saved to ~/%s'%fn)
        except Exception as e:
            print('Error: %s'%e)               
    
w_plot.on_click(on_plot_button_clicked)

def on_export_ass_button_clicked(b):
    try:
        background = collectionmosaic.select(0).add(15).divide(15)
        bgname = 'collectionmosaic'
        if w_collection.value == 'COPERNICUS/S1_GRD':  
            collection1 = getS2collection()
            count1 = collection1.size().getInfo()
            if count1>0:
        #      use sentinel-2 as video background if available                       
                background = ee.Image(collection1.first()) \
                                       .clip(poly) \
                                       .select('B8') 
                timestamp = background.get('system:time_start').getInfo()
                timestamp = time.gmtime(int(timestamp)/1000)
                timestamp = time.strftime('%x', timestamp)
                bgname = 'sentinel-2 '+ str(timestamp) 
                background = background.divide(5000)
        cmaps = ee.Image.cat(cmap,smap,fmap,bmap,background).rename(['cmap','smap','fmap']+timestamplist1[1:]+['background'])  
        assexport = ee.batch.Export.image.toAsset(cmaps.clip(poly),
                                    description='assetExportTask', 
                                    assetId=w_exportassetsname.value,scale=w_exportscale.value,maxPixels=1e9)      
        assexport.start()
        with w_out: 
            w_out.clear_output() 
            print('Exporting change maps to %s\n task id: %s'%(w_exportassetsname.value,str(assexport.id)))
    #  export metadata to drive
        if w_collection.value == 'COPERNICUS/S1_GRD': 
            times = [timestamp[1:9] for timestamp in timestamplist1]
            metadata = ee.List(['SEQUENTIAL OMNIBUS: '+time.asctime(),  
                                'Collection: '+w_collection.value,
                                'Asset export name: '+w_exportassetsname.value,  
                                'ENL: '+str(w_enl.value),
                                'Export scale (m): '+str(w_exportscale.value),
                                'Nominal scale (m): '+str(cmap.projection().nominalScale().getInfo()),
                                'Orbit pass: '+w_orbitpass.value,    
                                'Significance: '+str(w_significance.value),  
                                'Series length: '+str(len(times)),
                                'Timestamps: '+str(times)[1:-1],
                                'Rel orbit numbers: '+str(rons),
                                'Platform: '+w_platform.value,
                                'Background image: '+bgname,
                                'Mean incidence angles: '+str(mean_incidence),
                                'Used 5x5 median filter: '+str(w_median.value)]) \
                                .cat(['Polygon:']) \
                                .cat(poly.getInfo()['coordinates'][0]) 
        else:
            metadata = ee.List(['SEQUENTIAL OMNIBUS: '+time.asctime(),  
                                'Collection: '+w_collection.value,
                                'Asset export name: '+w_exportassetsname.value,  
                                'ENL: '+str(w_enl.value),  
                                'Export scale (m): '+str(w_exportscale.value),
                                'Nominal scale (m): '+str(cmap.projection().nominalScale().getInfo()),
                                'Significance: '+str(w_significance.value),  
                                'Series length: '+str(count),
                                'Used 3x3 median filter: '+str(w_median.value)]) 
        fileNamePrefix=w_exportassetsname.value.replace('/','-')  
        gdexport = ee.batch.Export.table.toDrive(ee.FeatureCollection(metadata.map(makefeature)),
                             description='driveExportTask_meta', 
                             folder = 'gee',
                             fileNamePrefix=fileNamePrefix )        
        gdexport.start()
        with w_out:
            print('Exporting metadata to Drive/gee/%s\n task id: %s'%(fileNamePrefix,str(gdexport.id)))    
    except Exception as e:
        with w_out:
            print('Error: %s'%e)                                          
    
w_export_ass.on_click(on_export_ass_button_clicked) 

def on_export_drv_button_clicked(b):
    try:
        cmaps = ee.Image.cat(cmap,smap,fmap,bmap).rename(['cmap','smap','fmap']+timestamplist1[1:])  
        fileNamePrefix=w_drive.value.replace('/','-')            
        gdexport = ee.batch.Export.image.toDrive(cmaps.byte().clip(poly),
                                    description='driveExportTask', 
                                    folder = 'gee',
                                    fileNamePrefix=fileNamePrefix,scale=w_exportscale.value,maxPixels=1e9)   
        gdexport.start()
        with w_out:
            w_out.clear_output()
            print('Exporting change maps to Drive/gee/%s\n task id: %s'%(fileNamePrefix,str(gdexport.id))) 

#      export metadata to drive
        if w_collection.value == 'COPERNICUS/S1_GRD': 
            times = [timestamp[1:9] for timestamp in timestamplist1]
            metadata = ee.List(['SEQUENTIAL OMNIBUS: '+time.asctime(),  
                                'Collection: '+w_collection.value,
                                'Drive export name: '+w_drive.value,  
                                'ENL: '+str(w_enl.value),
                                'Export scale (m): '+str(w_exportscale.value),
                                'Nominal scale (m): '+str(cmap.projection().nominalScale().getInfo()),
                                'Orbit pass: '+w_orbitpass.value,    
                                'Significance: '+str(w_significance.value),  
                                'Series length: '+str(len(times)),
                                'Timestamps: '+str(times)[1:-1],
                                'Rel orbit number(s): '+str(rons),
                                'Platform: '+w_platform.value,
                                'Mean incidence angles: '+str(mean_incidence),
                                'Used 5x5 median filter: '+str(w_median.value)]) \
                                .cat(['Polygon:']) \
                                .cat(poly.getInfo()['coordinates'][0]) 
        else:
            metadata = ee.List(['SEQUENTIAL OMNIBUS: '+time.asctime(),  
                                'Collection: '+w_collection.value,
                                'Drive export name: '+w_drive.value,  
                                'ENL: '+str(w_enl.value),  
                                'Export scale (m): '+str(w_exportscale.value),
                                'Nominal scale (m): '+str(cmap.projection().nominalScale().getInfo()),
                                'Significance: '+str(w_significance.value),  
                                'Series length: '+str(count),
                                'Used 5x5 median filter: '+str(w_median.value)]) 
        fileNamePrefix=w_exportdrivename.value.replace('/','-')  
        gdexport = ee.batch.Export.table.toDrive(ee.FeatureCollection(metadata.map(makefeature)),
                             description='driveExportTask_meta', 
                             folder = 'gee',
                             fileNamePrefix=fileNamePrefix )
        gdexport.start()
        with w_out:
            print('Exporting metadata to Drive/gee/%s\n task id: %s'%(fileNamePrefix,str(gdexport.id)))                   
    except Exception as e:
        with w_out:
            print('Error: %s'%e) 

w_export_drv.on_click(on_export_drv_button_clicked) 
        
def on_export_atsf_button_clicked(b):
#  export last ATSF image to drive together with unfiltered, hybrid filtered version and log image   
    try:          
        img_last = ee.Image(ee.List(imList).get(-1))
        img_rl = ee.Image(refinedLee(img_last))   
        if w_maskwater.value:
            img_atsf = ee.Image(avimg).updateMask(watermask)
        else:
            img_atsf = ee.Image(avimg)  
        img_log = ee.Image(avimglog)                   
        img_hybrid = img_atsf.where(img_log.lt(ee.Number(count).divide(3)),img_rl)   
        with w_out:       
            w_out.clear_output()     
            print('Exporting ATSF (adaptive temporal speckle filter) image to Drive')            
            gdexport1 = ee.batch.Export.image.toDrive(img_atsf,
                                                      description='driveExportTask_atsf', 
                                                      folder = 'gee',
                                                      fileNamePrefix = timestamplist1[-1]+'_atsf',
                                                      crs = archive_crs,
                                                      scale = w_exportscale.value,
                                                      maxPixels = 1e10)
            gdexport1.start()    
            print('Exporting ATSF log image to Drive')
            gdexport2 = ee.batch.Export.image.toDrive(ee.Image(img_log),
                                                      description='driveExportTask_atsf_log', 
                                                      folder = 'gee',
                                                      fileNamePrefix = timestamplist1[-1]+'_atsf_log',
                                                      crs = archive_crs,
                                                      scale = w_exportscale.value,
                                                      maxPixels = 1e10)
            gdexport2.start()  
            if w_collection.value == 'COPERNICUS/S1_GRD':
                print('Exporting hybrid image to Drive')
                gdexport3 = ee.batch.Export.image.toDrive(ee.Image(img_hybrid),
                                                          description='driveExportTask_atsf_hybrid', 
                                                          folder = 'gee',
                                                          fileNamePrefix = timestamplist1[-1]+'_atsf_hybrid',
                                                          crs = archive_crs,
                                                          scale = w_exportscale.value,
                                                          maxPixels = 1e10)
                gdexport3.start()  
            print('Exporting last image to Drive')
            gdexport4 = ee.batch.Export.image.toDrive(ee.Image(img_last),
                                                      description='driveExportTask_last', 
                                                      folder = 'gee',
                                                      fileNamePrefix = timestamplist1[-1],
                                                      crs = archive_crs,
                                                      scale = w_exportscale.value,
                                                      maxPixels = 1e10)
            gdexport4.start()   
# #          additionally export every kth ATSF filtered image to drive            
#             print('Exporting ATSF series to Drive ...')
#             k = 5
#             for i in range(0,count,k):
#                 if i<10:
#                     pad = '0'
#                 else:
#                     pad = ''
#                 image = ee.Image(avimgs.get(i)).clip(poly)
#                 gdexport5 = ee.batch.Export.image.toDrive(image,
#                                                           description='driveExportTask_series_'+pad+str(i), 
#                                                           folder = 'gee',
#                                                           fileNamePrefix = timestamplist1[i],
#                                                           crs = archive_crs,
#                                                           scale = w_exportscale.value,
#                                                           maxPixels = 1e10)     
#                gdexport5.start()                       
                                
    except Exception as e:
        with w_out:
            print('Error: %s'%e)        

w_export_atsf.on_click(on_export_atsf_button_clicked)       

def on_export_s2_button_clicked(b):
#  export clipped s2 image   
    try:           
        with w_out:       
            w_out.clear_output()     
            print('Exporting s2 image (optical bands B2, B3, B4 only) to Drive')            
            gdexport = ee.batch.Export.image.toDrive(s2_image.clip(poly),
                                                      description='driveExportTask_s2', 
                                                      folder = 'gee',
                                                      fileNamePrefix = 's2_%s_optical'%timestamps2,
                                                      crs = archive_crs,
                                                      scale = w_exportscale.value,
                                                      maxPixels = 1e11)
            gdexport.start()
    except Exception as e:
        with w_out:
            print('Error: %s'%e)        
            
w_export_s2.on_click(on_export_s2_button_clicked)                     
                          
def run():
    global m,center
#    center = list(reversed(poly.centroid().coordinates().getInfo()))
    center = [51.0,6.4]
    osm = basemap_to_tiles(basemaps.OpenStreetMap.Mapnik)
    ews = basemap_to_tiles(basemaps.Esri.WorldStreetMap)
    ewi = basemap_to_tiles(basemaps.Esri.WorldImagery)
    
    dc = DrawControl(polyline={},circlemarker={})
    dc.rectangle = {"shapeOptions": {"fillColor": "#0000ff","color": "#0000ff","fillOpacity": 0.05}}
    dc.polygon = {"shapeOptions": {"fillColor": "#0000ff","color": "#0000ff","fillOpacity": 0.05}}
    dc.on_draw(handle_draw)
    
    lc = LayersControl(position='topright')
    fs = FullScreenControl(position='topleft')
    mc = MeasureControl(position='topright',primary_length_unit = 'kilometers')

    m = Map(center=center, zoom=11, layout={'height':'500px'},layers=(ewi,ews,osm),controls=(mc,dc,lc,fs))   
#    m = Map(center=center, zoom=11, layout={'height':'500px'},controls=(lc,dc,fs,mc,sm_control)) 

    with w_out:
        w_out.clear_output()
        print('Algorithm output')
    display(m) 
    return box
    
