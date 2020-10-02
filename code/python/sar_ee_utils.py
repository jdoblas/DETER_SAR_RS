# -*- coding: utf-8 -*-
"""
collection of tools to obtain and process
SAR data using Google Earth Engine
@author: Juan Doblas (jdoblas@inpe.br)
"""
import os,sys,ee,math,time,datetime
def getS1dataFloat(date1,date2):
  return ee.ImageCollection("COPERNICUS/S1_GRD_FLOAT") \
  .filterMetadata('instrumentMode','equals','IW') \
  .filterMetadata('orbitProperties_pass','equals','DESCENDING') \
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')) \
  .select('VV','VH','angle') \
  .filterDate(date1,date2)
  
def getS1dataFloatByIngestionDate(date1,date2):
    date1p=date1.millis().multiply(1000)
    date2p=date2.millis().multiply(1000)
    return ee.ImageCollection("COPERNICUS/S1_GRD_FLOAT") \
    .filterMetadata('instrumentMode','equals','IW') \
    .filterMetadata('orbitProperties_pass','equals','DESCENDING') \
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')) \
    .select('VV','VH','angle') \
    .filterMetadata('system:version','greater_than',date1p)\
    .filterMetadata('system:version','not_greater_than',date2p)
    
def compute_pol_area(pol,proj='EPSG:5880'):
    area_ha=pol.area(maxError=1,proj=ee.Projection(proj)).divide(10000)
    return pol.set('area_ha',area_ha)

def execTask(task,check_interval=15,exit_on_error= True):
    """
    This function will execute and time a earth engine task.
    It is meant to stop execution of the main script until task ends
    
    Args:
        task (ee.task): The task to be executed
        check_interval (INT, optional): Interval to check task status. Defaults to 15.
        exit_on_error (BOOLEAN, optional): Wether script should terminate if task fails . Defaults to True.

    Returns:
        None.

    """
    
    start_time=time.time()
    task.start()
    task_id=task.status().get('id')
    task_status=task.status().get('state')
    start_time_f=datetime.datetime.now().strftime("%b %d %Y %H:%M:%S")
    print (f"Task {task_id} started at {start_time_f}")
    print ("waiting for task to complete, checking every "+str(check_interval)+" seconds")
    while task_status!='COMPLETED':
        #print(".",end="\r")
        time.sleep(check_interval)
        task_status=task.status().get('state')
        if (task_status=="FAILED" or task_status=='CANCELED'):
            print ()
            print ("Sorry, error on task:")
            print (task.status().get('error_message'))
            if exit_on_error:
                sys.exit()
    end_time=time.time()
    end_time_f=datetime.datetime.now().strftime("%b %d %Y %H:%M:%S")
    print ('')
    print (f"Task finished at {end_time_f}")
    print ("Duration: "+str(round((end_time-start_time)/60,2))+" minutes")

def toNatural(img):
  return ee.Image(ee.Image(10.0).pow(img.divide(10.0)).copyProperties(img,['system:time_start','sliceNumber']))

def toDB(img):
  return ee.Image(ee.Image(img).log10().multiply(10.0).copyProperties(img,['system:time_start','sliceNumber']))

def extractDates(col):
  datesList=ee.List(col.aggregate_array('system:time_start'))
  datesListFmt=datesList.map(formatDate)
  return datesListFmt

def formatDate(date):
  return ee.Date(date).format('YYYY-MM-dd')
  

def getNormalDistPdf(img,params):
  mean=params.select(0)
  sd=params.select(1)
  a=img.subtract(mean).divide(sd).pow(2).multiply(-0.5)
  return a.exp().multiply(0.39894228).divide(sd)
    
# Convert sigma0 db image to Gamma nought
def toGamma0(img):
  # WARNING: INPUT VALUES MUST BE IN DB!!!! BAND0: VV BAND1:VH
  img=img.addBands(getLIA(img).rename('LIA')) # Computes Local Incidence Angle (LIA) band
  lia=img.select('LIA')
  vv_gamma0=img.select(0).subtract(lia.multiply(math.pi/180.0).cos().log10().multiply(10.0))
  vh_gamma0=img.select(1).subtract(lia.multiply(math.pi/180.0).cos().log10().multiply(10.0))
  return img.addBands(vv_gamma0.rename('VVg0')).addBands(vh_gamma0.rename('VHg0'))

def toGamma0natural(img):
  # WARNING: INPUT VALUES MUST BE IN linear scale!!!! BAND0: VV BAND1:VH
  img=img.addBands(getLIA(img).rename('LIA')) # Computes Local Incidence Angle (LIA) band
  lia=img.select('LIA')
  vv_gamma0=img.select(0).divide(lia.multiply(math.pi/180.0).cos())
  vh_gamma0=img.select(1).divide(lia.multiply(math.pi/180.0).cos())
  return img.addBands(vv_gamma0.rename('VVg0')).addBands(vh_gamma0.rename('VHg0')).addBands(lia.rename('LIA'))

#/ Compute local incidence angle    #/
#/ by Felix Greifeneder, Guido Lemoine #/
def getLIA(img):
  srtm=ee.Image('USGS/SRTMGL1_003') # Loads MDT
  s1_inc = img.select('angle')
  s1_azimuth = ee.Terrain.aspect(s1_inc) \
                             .reduceRegion(ee.Reducer.mean(), s1_inc.get('system:footprint'), 100) \
                             .get('aspect')
  azimuthEdge = getDESCCorners(img)
  TrueAzimuth = azimuthEdge.get('azimuth') ;  # This should be some degree off the South direction (180), due to Earth rotation
  rotationFromSouth = ee.Number(TrueAzimuth).subtract(180.0)
  s1_azimuth = ee.Number(s1_azimuth).add(rotationFromSouth)
  srtm_slope = ee.Terrain.slope(srtm).select('slope')
  srtm_aspect = ee.Terrain.aspect(srtm).select('aspect')
  slope_projected = srtm_slope.multiply(ee.Image.constant(TrueAzimuth).subtract(90.0).subtract(srtm_aspect).multiply(math.pi/180).cos())
  lia = s1_inc.subtract(ee.Image.constant(90).subtract(ee.Image.constant(90).subtract(slope_projected))).abs()
  return lia

# Calculate True azimuth direction for  the near range image edge
def getDESCCorners(f):
  # Get the coords as a transposed array
  coords = ee.Array(f.geometry().coordinates().get(0)).transpose()
  crdLons = ee.List(coords.toList().get(0))
  crdLats = ee.List(coords.toList().get(1))
  minLon = crdLons.sort().get(0)
  maxLon = crdLons.sort().get(-1)
  minLat = crdLats.sort().get(0)
  maxLat = crdLats.sort().get(-1)
  azimuth = ee.Number(crdLons.get(crdLats.indexOf(minLat))).subtract(minLon) \
  .atan2(ee.Number(crdLats.get(crdLons.indexOf(minLon))).subtract(minLat)) \
  .multiply(180.0/math.pi) \
  .add(180.0);
  return ee.Feature(ee.Geometry.LineString([crdLons.get(crdLats.indexOf(maxLat)), maxLat,
    minLon, crdLats.get(crdLons.indexOf(minLon))]), { 'azimuth': azimuth}).copyProperties(f)

def makeFrostFilter(frostDamp, ksize):
    r=(ksize-1)/2
    # Frost filter by Guido Lemoine
    # Make sure to pass the frost damp factor (frostDamp) as a NEGATIVE number
    def Ffilter(image):
        distance_kernel = ee.Kernel.euclidean(r)
        # Square kernel
        kernel = ee.Kernel.square(r, 'pixels', False)
        # Get mean and variance
        mean = image.reduceNeighborhood(ee.Reducer.mean(), kernel)
        variance = image.reduceNeighborhood(ee.Reducer.variance(), kernel)
        B = variance.multiply(frostDamp).divide(mean.multiply(mean))
        W = B.exp().reduceNeighborhood(ee.Reducer.mean(), distance_kernel)
        return ee.Image(image.multiply(W).reduceNeighborhood(ee.Reducer.sum(), kernel).divide(W.reduceNeighborhood(ee.Reducer.sum(), kernel)).copyProperties(image,['system:time_start']))
    return Ffilter


def QueganYuFilter(imgCol,func):
  # This function will return a filtered collection
  # the filter is the proposed by Quegan&Yu (2001),
  # and uses a custom filter to find <J>
  # Warning: will only work on natural numbers (no DB).
  imgColMedian = imgCol.map(func)
  correctionFactorCol=imgCol.map(lambda img: img.divide(func(img)))
  correctionFactor=correctionFactorCol.sum().divide(imgCol.count())
  return imgColMedian.map(lambda img: img.multiply(correctionFactor).copyProperties(img,['system:time_start']))


def makeMedianFilter(ksize):
  def medianF(img):
    kernelMedian=ee.Kernel.square((ksize-1)/2,'pixels',False)
    return ee.Image(img.reduceNeighborhood(ee.Reducer.median(),kernelMedian)\
    .copyProperties(img,['system:time_start','sliceNumber']))
  return medianF

def refinedLeeFilter(img):
    # by Guido Lemoine
    # img must be in natural units, i.e. not in dB!
    img=ee.Image(img)
    # Set up 3x3 kernels 
    weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
    kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, False);
    
    mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
    variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3)
    # Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
    sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);
    sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, False);
    
    # Calculate mean and variance for the sampled windows and store as 9 bands
    sample_mean = mean3.neighborhoodToBands(sample_kernel); 
    sample_var = variance3.neighborhoodToBands(sample_kernel);
    
    # Determine the 4 gradients for the sampled windows
    gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
    gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
    gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
    gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());
    
    # And find the maximum gradient amongst gradient bands
    max_gradient = gradients.reduce(ee.Reducer.max());
    
    # Create a mask for band pixels that are the maximum gradient
    gradmask = gradients.eq(max_gradient);
    
    # duplicate gradmask bands: each gradient represents 2 directions
    gradmask = gradmask.addBands(gradmask);
    
    # Determine the 8 directions
    directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1)
    directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2))
    directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3))
    directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4))
    # The next 4 are the not() of the previous 4
    directions = directions.addBands(directions.select(0).Not().multiply(5));
    directions = directions.addBands(directions.select(1).Not().multiply(6));
    directions = directions.addBands(directions.select(2).Not().multiply(7));
    directions = directions.addBands(directions.select(3).Not().multiply(8));
    #Mask all values that are not 1-8
    directions = directions.updateMask(gradmask);
    # "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
    directions = directions.reduce(ee.Reducer.sum());  
    
    #var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
    #Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);
    
    sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));
    
    # Calculate localNoiseVariance
    sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);
    
    # Set up the 7*7 kernels for directional statistics
    rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));
    
    diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
      [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);
    
    rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, False);
    diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, False);
    
    # Create stacks for mean and variance using the original kernels. Mask with relevant direction.
    dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
    dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));
    
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));
    
    # and add the bands for rotated kernels
    for i in [1,2,3]:
      dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
      dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
      dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
      dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    
    # "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
    dir_mean = dir_mean.reduce(ee.Reducer.sum());
    dir_var = dir_var.reduce(ee.Reducer.sum());
    
    # A finally generate the filtered value
    varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));
    
    b = varX.divide(dir_var);
    
    result = dir_mean.add(b.multiply(img.subtract(dir_mean)));
    return result.arrayFlatten([['sum']]).rename(img.bandNames()).copyProperties(img,['system:time_start']);



def ee_export_vector_silent(ee_object, filename, selectors=None):
    """Silent version of ee_export_vector of Qiusheng Wu
    Exports Earth Engine FeatureCollection to other formats, including shp, csv, json, kml, and kmz.
    Args:
        ee_object (object): ee.FeatureCollection to export.
        filename (str): Output file name.
        selectors (list, optional): A list of attributes to export. Defaults to None.
    """
    import requests
    import zipfile
    #ee_initialize()

    if not isinstance(ee_object, ee.FeatureCollection):
        print('The ee_object must be an ee.FeatureCollection.')
        return

    allowed_formats = ['csv', 'json', 'kml', 'kmz', 'shp']
    filename = os.path.abspath(filename)
    basename = os.path.basename(filename)
    name = os.path.splitext(basename)[0]
    filetype = os.path.splitext(basename)[1][1:].lower()
    filename_shp = filename

    if filetype == 'shp':
        filename = filename.replace('.shp', '.zip')

    if not (filetype.lower() in allowed_formats):
        print('The file type must be one of the following: {}'.format(
            ', '.join(allowed_formats)))
        return

    if selectors is None:
        selectors = ee_object.first().propertyNames().getInfo()
    elif not isinstance(selectors, list):
        print("selectors must be a list, such as ['attribute1', 'attribute2']")
        return
    else:
        allowed_attributes = ee_object.first().propertyNames().getInfo()
        for attribute in selectors:
            if not (attribute in allowed_attributes):
                print('Attributes must be one chosen from: {} '.format(
                    ', '.join(allowed_attributes)))
                return

    try:
        #print('Generating URL ...')
        url = ee_object.getDownloadURL(
            filetype=filetype, selectors=selectors, filename=name)
        #print('Downloading data from {}\nPlease wait ...'.format(url))
        r = requests.get(url, stream=True)

        if r.status_code != 200:
            print('An error occurred while downloading. \n Retrying ...')
            try:
                new_ee_object = ee_object#.map(filter_polygons)
                print('Generating URL ...')
                url = new_ee_object.getDownloadURL(
                    filetype=filetype, selectors=selectors, filename=name)
                print('Downloading data from {}\nPlease wait ...'.format(url))
                r = requests.get(url, stream=True)
            except Exception as e:
                print(e)

        with open(filename, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=1024):
                fd.write(chunk)
    except Exception as e:
        print('An error occurred while downloading.')
        print(e)
        return

# =============================================================================
#     try:
#         if filetype == 'shp':
#             z = zipfile.ZipFile(filename)
#             z.extractall(os.path.dirname(filename))
#             #os.remove(filename)
#             #filename = filename.replace('.zip', '.shp')
# 
#         #print('Data downloaded to {}'.format(filename))
#     except Exception as e:
#         print(e)
# 
# =============================================================================
