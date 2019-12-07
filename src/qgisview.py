# **************************************************
# Viewer for Exported Sequential Omnibus Change Maps
# **************************************************

import ee
from ee_plugin import Map

#util = require('users/mortcanty/changedetection:utilities');
#text = require('users/gena/packages:text');

image = ee.Image('users/mortcanty/omnibus/buziflood')

k = image.bandNames().length().subtract(4).getInfo();
k1 = k+3;
jet = ['black','blue','cyan', 'yellow','red'];
rgy = ['black','red','green','yellow'] 
vis = {'min':0, 'max':k, 'palette':jet};
vis1 = {'min':0, 'max':3, 'palette':rgy};

Map.centerObject(image,13);
#Map.add(util.makeLegend(vis));

Map.addLayer(image.select('fmap').multiply(2),vis,'fmap*2');
Map.addLayer(image.select('smap'),vis,'smap');
Map.addLayer(image.select('cmap'),vis,'cmap');
Map.addLayer(image.select('background'),{'min':0,'max':1},'background');