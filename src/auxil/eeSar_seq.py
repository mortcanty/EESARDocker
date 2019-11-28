'''
Created on 21.06.2018

@author: mort

ipywidget interface to the GEE for sequential SAR change detection

'''
import ee, time, warnings, math
import ipywidgets as widgets
from IPython.display import display
from ipyleaflet import (Map,DrawControl,TileLayer,
                        basemaps,basemap_to_tiles,
                        LayersControl,
                        SplitMapControl)
from auxil.eeWishart import omnibus

ee.Initialize()

warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

poly = ee.Geometry.Polygon([[6.30154, 50.948329], [6.293307, 50.877329], 
                            [6.427091, 50.875595], [6.417486, 50.947464], 
                            [6.30154, 50.948329]])

center = list(reversed(poly.centroid().coordinates().getInfo()))

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
    ''' clip a list of images '''
    imlist = ee.List(ee.Dictionary(prev).get('imlist'))
    poly = ee.Dictionary(prev).get('poly') 
    enl = ee.Number(ee.Dictionary(prev).get('enl'))   
    imlist = imlist.add(ee.Image(current).multiply(enl).clip(poly))
    return ee.Dictionary({'imlist':imlist,'poly':poly,'enl':enl})

def makefeature(data):
    ''' for exporting as CSV to Drive '''
    return ee.Feature(None, {'data': data})

def handle_draw(self, action, geo_json):
    global poly
    if action == 'created':
        coords =  geo_json['geometry']['coordinates']
        poly = ee.Geometry.Polygon(coords)
        w_preview.disabled = True
        w_export_ass.disabled = True
        
dc = DrawControl(polyline={},circlemarker={})
dc.rectangle = {"shapeOptions": {"fillColor": "#0000ff","color": "#0000ff","fillOpacity": 0.1}}
dc.on_draw(handle_draw)

def GetTileLayerUrl(ee_image_object):
    map_id = ee.Image(ee_image_object).getMapId()
    tile_url_template = "https://earthengine.googleapis.com/map/{mapid}/{{z}}/{{x}}/{{y}}?token={token}"
    return tile_url_template.format(**map_id)

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
w_bmap = widgets.IntText(
    layout = widgets.Layout(width='50px'),
    value='1',
    description='',
    disabled=False
)
w_platform = widgets.RadioButtons(
    layout = widgets.Layout(width='20%'),
    options=['Both','A','B'],
     value='Both',
    description='Platform:',
    disabled=False
)
w_relativeorbitnumber = widgets.IntText(
    value='0',
    description='Rel orbit:',
    disabled=False
)
w_exportassetsname = widgets.Text(
    value='users/<username>/<path>',
    placeholder=' ',
    disabled=False
)
w_exportdrivename = widgets.Text(
    value='EarthEngineImages/<path>',
    placeholder=' ',
    disabled=False
)
w_exportscale = widgets.Text(
    layout = widgets.Layout(width='150px'),
    value='10',
    placeholder=' ',
    description='Scale ',
    disabled=False
)
w_series = widgets.Checkbox(
    value=False,
    description='Export series to drive',
    disabled=False
)
w_startdate = widgets.Text(
    value='2018-04-01',
    placeholder=' ',
    description='Start date:',
    disabled=False
)
w_enddate = widgets.Text(
    value='2018-10-01',
    placeholder=' ',
    description='End date:',
    disabled=False
)
w_median = widgets.Checkbox(
    value=True,
    description='3x3 Median filter',
    disabled=False
)
w_S2 = widgets.Checkbox(
    value=True,
    description='Show best S2',
    disabled=False
)
w_Q = widgets.Checkbox(
    value=True,
    description='Quick Preview',
    disabled=False
)
w_significance = widgets.BoundedFloatText(
    layout = widgets.Layout(width='200px'),
    value='0.0001',
    min=0,
    max=0.05,
    step=0.0001,
    description='Significance:',
    disabled=False
)
w_opacity = widgets.BoundedFloatText(
    layout = widgets.Layout(width='50%'),
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
w_maskland = widgets.Checkbox(
    value=True,
    description='Water mask',
    disabled=False
)
w_text = widgets.Textarea(
    layout = widgets.Layout(width='99%'),
    value = 'Algorithm output',
    rows = 4,
    disabled = False
)

w_collect = widgets.HBox([w_collection,w_enl])
w_opac = widgets.VBox([w_opacity,w_maskchange,w_maskland],layout = widgets.Layout(width='50%'))
w_run = widgets.Button(description="Collect")
w_preview = widgets.Button(description="Preview",disabled=True)
w_export_ass = widgets.Button(description='ExportToAssets',disabled=True)
w_export_drv = widgets.Button(description='ExportToDrive',disabled=True)
w_dates = widgets.HBox([w_relativeorbitnumber,w_startdate,w_enddate])
w_change = widgets.HBox([w_changemap,w_bmap])
w_orbit = widgets.HBox([w_orbitpass,w_platform,w_change,w_opac])
w_exp = widgets.HBox([w_export_ass,w_exportassetsname,w_export_drv,w_exportdrivename,w_series])
w_signif = widgets.HBox([w_significance,w_S2,w_Q,w_median,w_exportscale],layout = widgets.Layout(width='99%'))
w_rse = widgets.HBox([w_run,w_preview])

box = widgets.VBox([w_text,w_collect,w_dates,w_orbit,w_signif,w_rse,w_exp])

def on_widget_change(b):
    w_preview.disabled = True
    w_export_ass.disabled = True
    w_export_drv.disabled = True

w_orbitpass.observe(on_widget_change,names='value')
w_platform.observe(on_widget_change,names='value')
w_relativeorbitnumber.observe(on_widget_change,names='value')
w_startdate.observe(on_widget_change,names='value')
w_collection.observe(on_widget_change,names='value')
w_enl.observe(on_widget_change,names='value')
w_enddate.observe(on_widget_change,names='value')
w_median.observe(on_widget_change,names='value')
w_significance.observe(on_widget_change,names='value')

def getS1collection(coords):
    return ee.ImageCollection('COPERNICUS/S1_GRD') \
                      .filterBounds(ee.Geometry.Point(coords.get(0))) \
                      .filterBounds(ee.Geometry.Point(coords.get(1))) \
                      .filterBounds(ee.Geometry.Point(coords.get(2))) \
                      .filterBounds(ee.Geometry.Point(coords.get(3))) \
                      .filterDate(ee.Date(w_startdate.value), ee.Date(w_enddate.value)) \
                      .filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV','VH'])) \
                      .filter(ee.Filter.eq('resolution_meters', 10)) \
                      .filter(ee.Filter.eq('instrumentMode', 'IW')) \
                      .filter(ee.Filter.eq('orbitProperties_pass', w_orbitpass.value))     
                      
def getS2collection(coords):
    return ee.ImageCollection('COPERNICUS/S2') \
                      .filterBounds(ee.Geometry.Point(coords.get(0))) \
                      .filterBounds(ee.Geometry.Point(coords.get(1))) \
                      .filterBounds(ee.Geometry.Point(coords.get(2))) \
                      .filterBounds(ee.Geometry.Point(coords.get(3))) \
                      .filterDate(ee.Date(w_startdate.value),ee.Date(w_enddate.value)) \
                      .sort('CLOUDY_PIXEL_PERCENTAGE',True)
 

def on_run_button_clicked(b):
    global result,m,collection,poly,imList,count,timestamplist1, \
           w_startdate,w_enddate,w_orbitpass,w_changemap,s2_image, \
           w_relativeorbitnumber,w_significance,w_median,w_Q, \
           mean_incidence,collectionmean,archive_crs,coords
    try:       
        if w_collection.value == 'COPERNICUS/S1_GRD':    
            w_text.value = 'running on GEE archive ...' 
            coords = ee.List(poly.bounds().coordinates().get(0))
            collection = getS1collection(coords) 
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
    #      make timestamps in YYYYMMDD format            
            timestamplist = [x.replace('/','') for x in timestamplist]  
            timestamplist = ['T20'+x[4:]+x[0:4] for x in timestamplist]
    #      in case of duplicates add running integer
            timestamplist1 = [timestamplist[i] + '_' + str(i+1) for i in range(len(timestamplist))]    
            relativeorbitnumbers = map(int,ee.List(collection.aggregate_array('relativeOrbitNumber_start')).getInfo())
            rons = list(set(relativeorbitnumbers))
            txt = 'running on GEE archive ...'
            txt += 'Images found: %i, platform: %s \n'%(count,w_platform.value)
            txt += 'Number of 10m pixels contained: %i \n'%math.floor(poly.area().getInfo()/100.0)
            txt += 'Acquisition dates: %s\n'%str(timestamplist)
            txt += 'Relative orbit numbers: '+str(rons)+'\n'
            if len(rons)==1:
                mean_incidence = get_incidence_angle(collection.first())
                txt += 'Mean incidence angle: %f'%mean_incidence
            else:
                txt += 'Mean incidence angle: (select one rel. orbit)'
            pcollection = collection.map(get_vvvh)
            pList = pcollection.toList(100)   
            first = ee.Dictionary({'imlist':ee.List([]),'poly':poly,'enl':ee.Number(w_enl.value)}) 
            imList = ee.Dictionary(pList.iterate(clipList,first)).get('imlist')
            collectionmean = collection.mean().select(0).clip(poly).rename('b0')     
            percentiles = collectionmean.reduceRegion(ee.Reducer.percentile([2,98]),scale=float(w_exportscale.value),maxPixels=10e9)
            mn = ee.Number(percentiles.get('b0_p2'))
            mx = ee.Number(percentiles.get('b0_p98'))        
            vorschau = collectionmean.visualize(min=mn, max=mx, opacity=w_opacity.value) 
        else:
            txt = 'running on local collection %s ...\n'%w_collection.value 
            collection = ee.ImageCollection(w_collection.value)
            count = collection.size().getInfo()  
            txt += 'Images found: %i\n'%count           
            collectionfirst = ee.Image(collection.first())  
            poly = collectionfirst.geometry()   
            coords = ee.List(poly.bounds().coordinates().get(0))   
            center = poly.centroid().coordinates().getInfo()
            center.reverse()
            m.center = center                
            w_exportscale.value = str(collectionfirst.projection().nominalScale().getInfo()) 
            if collectionfirst.get('system:time_start').getInfo() is not None:
                acquisition_times = ee.List(collection.aggregate_array('system:time_start')).getInfo()  
                timestamplist1 = []
                for timestamp in acquisition_times:
                    tmp = time.gmtime(int(timestamp)/1000)
                    timestamplist1.append(time.strftime('%x', tmp))            
                timestamplist1 = [x.replace('/','') for x in timestamplist1]  
                timestamplist1 = ['T20'+x[4:]+x[0:4] for x in timestamplist1]    
                txt += 'Acquisition dates: %s'%str(timestamplist1)    
            else:
                timestamplist1 = ['T%i'%(i+1) for i in range(count)]
                txt += 'No time property available: acquisitions: %s'%str(timestamplist1)
            collectionmean = collection.mean().select(0)
            percentiles = collectionfirst.select(0).rename('b0').reduceRegion(ee.Reducer.percentile([0,98]),maxPixels=10e9)
            mn = ee.Number(percentiles.get('b0_p0'))
            mx = ee.Number(percentiles.get('b0_p98'))        
            vorschau = collectionmean.visualize(min=mn, max=mx, opacity=w_opacity.value)           
            imList = collection.toList(100)
#      get GEE S1 archive crs for eventual image series export               
        archive_crs = ee.Image(getS1collection(coords).first()).select(0).projection().crs().getInfo()
#      run the algorithm        
        result = ee.Dictionary(omnibus(imList,w_significance.value,w_enl.value,w_median.value))
        w_preview.disabled = False
        s2_image = None
#      display collection or S2                
        if len(m.layers)>1:
            m.remove_layer(m.layers[1])
        if w_S2.value:
#          display sentinel-2 if available              
            collection2 = getS2collection(coords) 
            count1 = collection2.size().getInfo()
            if count1>0:    
                s2_image =  ee.Image(collection2.first()).select(['B8','B4','B3']).clip(poly)        
                percentiles = s2_image.reduceRegion(ee.Reducer.percentile([2,98]),scale=float(w_exportscale.value),maxPixels=10e9)         
                mn = percentiles.values(['B8_p2','B4_p2','B3_p2'])
                mx = percentiles.values(['B8_p98','B4_p98','B3_p98'])
                vorschau = s2_image.visualize(min=mn,max=mx,opacity=w_opacity.value)           
                timestamp = s2_image.get('system:time_start').getInfo() 
                timestamp = time.gmtime(int(timestamp)/1000)
                timestamp = time.strftime('%x', timestamp).replace('/','')
                timestamp = '20'+timestamp[4:]+timestamp[0:4]
                txt += '\nSentinel-2 from %s'%timestamp
        w_text.value = txt           
        m.add_layer(TileLayer(url=GetTileLayerUrl(vorschau)))
    except Exception as e:
        w_text.value =  'Error: %s'%e

w_run.on_click(on_run_button_clicked)

def on_preview_button_clicked(b):
    global result,m,cmap,smap,fmap,bmap,w_changemap
    
    landmask = ee.Image('UMD/hansen/global_forest_change_2015').select('datamask').eq(1)    
    
    jet = 'black,blue,cyan,yellow,red'
    rgy = 'black,red,green,yellow'
    smap = ee.Image(result.get('smap')).byte()
    cmap = ee.Image(result.get('cmap')).byte()
    fmap = ee.Image(result.get('fmap')).byte() 
    bmap = ee.Image(result.get('bmap')).byte()
    
    cmaps = ee.Image.cat(cmap,smap,fmap,bmap).rename(['cmap','smap','fmap']+timestamplist1[1:])
    downloadpath = cmaps.getDownloadUrl({'scale':w_exportscale.value}) 
    txt = 'Download change maps from this URL (may be unreliable):\n'
    txt += downloadpath+'\n'
    palette = jet
    if w_changemap.value=='First':
        mp = smap
        mx = count
        txt += 'Interval of first change:\n blue = early, red = late'
    elif w_changemap.value=='Last':
        mp=cmap
        mx = count
        txt += 'Interval of last change:\n blue = early, red = late'
    elif w_changemap.value=='Frequency':
        mp = fmap
        mx = count/2
        txt += 'Interval of change frequency :\n blue = few, red = many'
    else:
        sel = int(w_bmap.value)
        sel = min(sel,count-1)
        sel = max(sel,1)
        txt = 'Bitemporal: %s-->%s\n'%(timestamplist1[sel-1],timestamplist1[sel])
        txt += 'red = positive definite, green = negative definite, yellow = indefinite'     
        mp = bmap.select(sel-1).clip(poly)
        palette = rgy
        mx = 3
    w_text.value = txt
    if len(m.layers)>1:
        m.remove_layer(m.layers[1])
    if not w_Q.value:
        mp = mp.reproject(crs=archive_crs,scale=float(w_exportscale.value))
    if w_maskland.value==True:
        mp = mp.updateMask(landmask)
    if w_maskchange.value==True:    
        mp = mp.updateMask(mp.gt(0))
    m.add_layer(TileLayer(url=GetTileLayerUrl(mp.visualize(min=0, max=mx, palette=palette,opacity = w_opacity.value))))
    w_export_ass.disabled = False
    w_export_drv.disabled = False
    
w_preview.on_click(on_preview_button_clicked)   

def on_export_ass_button_clicked(b):
    global w_exportassetsname
    background = collectionmean.add(15).divide(15)
    bgname = 'collectionmean'
    if w_collection.value == 'COPERNICUS/S1_GRD':  
        collection1 = getS2collection(coords)
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
            background = background.divide(10000)
    cmaps = ee.Image.cat(cmap,smap,fmap,bmap,background).rename(['cmap','smap','fmap']+timestamplist1[1:]+['background'])                
    assexport = ee.batch.Export.image.toAsset(cmaps,
                                description='assetExportTask', 
                                assetId=w_exportassetsname.value,scale=w_exportscale.value,maxPixels=1e9)
    assexportid = str(assexport.id)
    w_text.value= 'Exporting change maps to %s\n task id: %s'%(w_exportassetsname.value,assexportid)
    assexport.start()  
    
#  export metadata to drive
    if w_collection.value == 'COPERNICUS/S1_GRD': 
        times = [timestamp[1:9] for timestamp in timestamplist1]
        metadata = ee.List(['SEQUENTIAL OMNIBUS: '+time.asctime(),  
                            'Collection: '+w_collection.value,
                            'Asset export name: '+w_exportassetsname.value,  
                            'ENL: '+str(w_enl.value),
                            'Export scale (m): '+w_exportscale.value,
                            'Nominal scale (m): '+str(cmap.projection().nominalScale().getInfo()),
                            'Orbit pass: '+w_orbitpass.value,    
                            'Significance: '+str(w_significance.value),  
                            'Series length: '+str(len(times)),
                            'Timestamps: '+str(times)[1:-1],
                            'Rel orbit numbers: '+str(w_relativeorbitnumber.value),
                            'Platform: '+w_platform.value,
                            'Background image: '+bgname,
                            'Mean incidence angles: '+str(mean_incidence),
                            'Used 3x3 median filter: '+str(w_median.value)]) \
                            .cat(['Polygon:']) \
                            .cat(poly.getInfo()['coordinates'][0]) 
    else:
        metadata = ee.List(['SEQUENTIAL OMNIBUS: '+time.asctime(),  
                            'Collection: '+w_collection.value,
                            'Asset export name: '+w_exportassetsname.value,  
                            'ENL: '+str(w_enl.value),  
                            'Export scale (m): '+w_exportscale.value,
                            'Nominal scale (m): '+str(cmap.projection().nominalScale().getInfo()),
                            'Significance: '+str(w_significance.value),  
                            'Series length: '+str(count),
                            'Used 3x3 median filter: '+str(w_median.value)]) 
    fileNamePrefix=w_exportassetsname.value.replace('/','-')  
    gdexport = ee.batch.Export.table.toDrive(ee.FeatureCollection(metadata.map(makefeature)),
                         description='driveExportTask_meta', 
                         folder = 'EarthEngineImages',
                         fileNamePrefix=fileNamePrefix )
    w_text.value += '\n Exporting metadata to Drive/EarthEngineImages/%s\n task id: %s'%(fileNamePrefix,str(gdexport.id))            
    gdexport.start()         
    if w_series.value == True:
        imlist = ee.List(imList)
        w_text.value += '\n Exporting time series of %i images to Drive'%count
        for i in range(count):
            if i<10:
                pad = '0'
            else:
                pad = ''
            image = ee.Image(imlist.get(i))
            gdexport1 = ee.batch.Export.image.toDrive(image,
                                                      description='driveExportTask_series_'+pad+str(i), 
                                                      folder = 'EarthEngineImages',
                                                      fileNamePrefix = fileNamePrefix+'-'+timestamplist1[i],
                                                     # fileNamePrefix = fileNamePrefix+'_'+pad+str(i),
                                                      crs = archive_crs,
                                                      scale = w_exportscale.value,
                                                      maxPixels = 1e10)
            gdexport1.start()  
        if s2_image is not None:
            w_text.value += '\n Exporting s2 image to Drive'
            gdexport2 = ee.batch.Export.image.toDrive(s2_image,
                                                      description='driveExportTask_s2', 
                                                      folder = 'EarthEngineImages',
                                                      fileNamePrefix = fileNamePrefix+'_s2',
                                                      crs = archive_crs,
                                                      scale = w_exportscale.value,
                                                      maxPixels = 1e10)
            gdexport2.start()                                                            
    
w_export_ass.on_click(on_export_ass_button_clicked) 

def on_export_drv_button_clicked(b):
    global w_exportdrivename
    background = collectionmean.add(15).divide(15)
    bgname = 'collectionmean'
    if w_collection.value == 'COPERNICUS/S1_GRD':  
        collection1 = getS2collection(coords)
        count1 = collection1.size().getInfo()
        if count1>0:
    #      use sentinel-2 as video background if available                       
            background = ee.Image(collection1.first()) \
                                   .clip(poly) \
                                   .select('B8') 
            timestamp = background.get('system:time_start').getInfo()
            timestamp = time.gmtime(int(timestamp)/5000)
            timestamp = time.strftime('%x', timestamp)
            bgname = 'sentinel-2 '+ str(timestamp) 
            background = background.divide(10000)
    cmaps = ee.Image.cat(cmap,smap,fmap,bmap).rename(['cmap','smap','fmap']+timestamplist1[1:])  
    fileNamePrefix=w_exportdrivename.value.replace('/','-')            
    gdexport = ee.batch.Export.image.toDrive(cmaps.byte(),
                                description='driveExportTask', 
                                folder = 'EarthEngineImages',
                                fileNamePrefix=fileNamePrefix,scale=w_exportscale.value,maxPixels=1e9)   
    w_text.value= 'Exporting change maps to Drive/EarthEngineImages/%s\n task id: %s'%(fileNamePrefix,str(gdexport.id)) 
    gdexport.start()  
    
    
#  export metadata to drive
    if w_collection.value == 'COPERNICUS/S1_GRD': 
        times = [timestamp[1:9] for timestamp in timestamplist1]
        metadata = ee.List(['SEQUENTIAL OMNIBUS: '+time.asctime(),  
                            'Collection: '+w_collection.value,
                            'Drive export name: '+w_exportdrivename.value,  
                            'ENL: '+str(w_enl.value),
                            'Export scale (m): '+w_exportscale.value,
                            'Nominal scale (m): '+str(cmap.projection().nominalScale().getInfo()),
                            'Orbit pass: '+w_orbitpass.value,    
                            'Significance: '+str(w_significance.value),  
                            'Series length: '+str(len(times)),
                            'Timestamps: '+str(times)[1:-1],
                            'Rel orbit numbers: '+str(w_relativeorbitnumber.value),
                            'Platform: '+w_platform.value,
                            'Background image: '+bgname,
                            'Mean incidence angles: '+str(mean_incidence),
                            'Used 3x3 median filter: '+str(w_median.value)]) \
                            .cat(['Polygon:']) \
                            .cat(poly.getInfo()['coordinates'][0]) 
    else:
        metadata = ee.List(['SEQUENTIAL OMNIBUS: '+time.asctime(),  
                            'Collection: '+w_collection.value,
                            'Drive export name: '+w_exportdrivename.value,  
                            'ENL: '+str(w_enl.value),  
                            'Export scale (m): '+w_exportscale.value,
                            'Nominal scale (m): '+str(cmap.projection().nominalScale().getInfo()),
                            'Significance: '+str(w_significance.value),  
                            'Series length: '+str(count),
                            'Used 3x3 median filter: '+str(w_median.value)]) 
    fileNamePrefix=w_exportdrivename.value.replace('/','-')  
    gdexport = ee.batch.Export.table.toDrive(ee.FeatureCollection(metadata.map(makefeature)),
                         description='driveExportTask_meta', 
                         folder = 'EarthEngineImages',
                         fileNamePrefix=fileNamePrefix )
    w_text.value += '\n Exporting metadata to Drive/EarthEngineImages/%s\n task id: %s'%(fileNamePrefix,str(gdexport.id))            
    gdexport.start()         
    if w_series.value == True:
        imlist = ee.List(imList)
        w_text.value += '\n Exporting time series of %i images to Drive'%count
        for i in range(count):
            if i<10:
                pad = '0'
            else:
                pad = ''
            image = ee.Image(imlist.get(i))
            gdexport1 = ee.batch.Export.image.toDrive(image,
                                                      description='driveExportTask_series_'+pad+str(i), 
                                                      folder = 'EarthEngineImages',
                                                      fileNamePrefix = fileNamePrefix+'_'+pad+str(i),
                                                      crs = archive_crs,
                                                      scale = w_exportscale.value,
                                                      maxPixels = 1e10)
            gdexport1.start()  
        if s2_image is not None:
            w_text.value += '\n Exporting s2 image to Drive'
            gdexport2 = ee.batch.Export.image.toDrive(s2_image,
                                                      description='driveExportTask_s2', 
                                                      folder = 'EarthEngineImages',
                                                      fileNamePrefix = fileNamePrefix+'_s2',
                                                      crs = archive_crs,
                                                      scale = w_exportscale.value,
                                                      maxPixels = 1e10)
            gdexport2.start()                                                            
    
w_export_drv.on_click(on_export_drv_button_clicked) 

def run():
    global m,dc,center
    center = list(reversed(poly.centroid().coordinates().getInfo()))
    m = Map(center=center, zoom=11, layout={'height':'400px'})
    osm = basemap_to_tiles(basemaps.OpenStreetMap.HOT)
    mb = TileLayer(url="https://api.mapbox.com/styles/v1/mapbox/satellite-streets-v9/tiles/256/{z}/{x}/{y}?access_token=pk.eyJ1IjoibWNhbnR5IiwiYSI6ImNpcjRsMmJxazAwM3hoeW05aDA1cmNkNzMifQ.d2UbIugbQFk2lnU8uHwCsQ")
    sm_control = SplitMapControl(left_layer=osm,right_layer=mb)
    m.add_control(dc)
    m.add_control(sm_control)
    m.add_control(LayersControl())
    display(m)
    display(box)
    
