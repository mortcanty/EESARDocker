import ee
import requests

def download(img,path='/home/mort/Downloads/',name='download',bands='VV',scale=10):
    crs = img.projection().crs().getInfo()
    img_url = ee.Image(img).getDownloadURL({'name':name,'bands':bands,'scale':scale})
    res = requests.get(img_url, stream=True)
#  The response is the zipped image which will contain download.tif
    outzip = path+'image.zip'
    if res.headers['Content-Type'] == "application/zip":
        handle = open(outzip, "wb")
        for chunk in res.iter_content(chunk_size=512):
            if chunk:  # filter out keep-alive new chunks
                handle.write(chunk)
        handle.close()
        return True
    else:
        print("Unexpected response content-type ", res.headers['Content-Type']) 
        return False   
         
if __name__ == '__main__':
    pass        

# ee.Initialize()
#  
# s2 = ee.ImageCollection('COPERNICUS/S2_SR').\
#   filterBounds(ee.Geometry.Point(7.1901, 51.479)).\
#   filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 10.0).\
#   sort('system:time_start', False)
#  
# img = ee.Image(s2.first()) # The Last will be the First
# img_url = img.getDownloadURL({'name': 'nice_name', 'bands': "B8,B4,B3", 'scale': 100})
#  
# print(img_url)
#  
# res = requests.get(img_url, stream=True)
#  
# # The response is the zipped image which will contain download.tif
# outzip = "image.zip"
#  
# if res.headers['Content-Type'] == "application/zip":
#     handle = open(outzip, "wb")
#     for chunk in res.iter_content(chunk_size=512+):
#         if chunk:  # filter out keep-alive new chunks
#             handle.write(chunk)
#     handle.close()
# else:
#     print("Unexpected response content-type ", res.headers['Content-Type'])