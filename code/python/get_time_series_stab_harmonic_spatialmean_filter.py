# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 17:13:46 2020

@author: juand
"""
import ee
import datetime
import pandas
import math
import logging
from sar_ee_utils import toDB,getS1dataFloat,toGamma0natural,refinedLeeFilter,makeMedianFilter,QueganYuFilter,getNormalDistPdf,ee_export_vector_silent
from get_forest_mask import get_forest_mask
from stabilize_SAR_time_series import stabilize_SAR_time_series,stabilize_SAR_time_series_upscale
from timeit import default_timer as timer

def toDB(img):
    return img.log10().multiply(10).copyProperties(img,['system:time_start']) 

def get_point_timeseries_2stab_filter(pol,ee_point,ee_date_detect,yearly_mask,optical_mask,sar_mask,sar_tmp_mask,\
                        detection_threshold=3,generalScale=20,forestBuffer=100,searchDistance=5000,\
                        detectPeriodMonths=3,learningPeriodYears=3,stabilize=1,\
                        stabilize_scale=0,shrink=0,include_S1B=0,\
                        adaptative_threshold=None):
    logging.getLogger('googleapicliet.discovery_cache').setLevel(logging.ERROR)
    AOI=ee_point.buffer(10000)
    # Define relevant dates
    date2=ee_date_detect
    date3=date2.advance(detectPeriodMonths,'months')# Begin of detection (inclusive)
    date1=date2.advance(-1*learningPeriodYears,'years') # Begin of training colection (inclusive)
    date1str=date1.format('YYYY-MM-dd').getInfo()
    date2str=date2.format('YYYY-MM-dd').getInfo()
    date3str=date3.format('YYYY-MM-dd').getInfo()
    # First status messages
    print ("Processing point")
    print ("Learning period start: "+date1str)
    print ("Detection start: "+date2str)
    print ("Detection end: "+date3str)
    # Define S1 colection from EE
    colS1=getS1dataFloat(date1,date3).filterBounds(AOI)\
          .sort('system:time_start')\
          .map(lambda img: img.updateMask(img.gt(.00002))) \
          .map(toGamma0natural) \
          .select(["VHg0","VVg0","LIA"])
    LIA=colS1.select(["LIA"]).mean().reduceRegion('mean',ee_point,generalScale).getInfo().get("LIA")
    print(f"Local Incidence Angle: {LIA}")
    if not include_S1B:
        colS1=colS1.filterMetadata("platform_number","equals","A")
    # Leave only the most represented orbit in every pixel
    # add relative orbit number
    def addOrbitBand(img):
        orb=ee.Number(img.get('relativeOrbitNumber_start'))
        return img.addBands(ee.Image(orb).rename('Orbit').toInt())
    colS1=colS1.map(addOrbitBand)
    orbit_mode=colS1.select('Orbit').reduce(ee.Reducer.mode())
    colS1=colS1.map(lambda img: img.updateMask(img.select('Orbit').eq(orbit_mode))).select(["VHg0","VVg0","LIA"])
    colSize=colS1.size().getInfo()    
    print ("# of collected images: "+str(colSize))
    if colSize>0:
        original_time_series=colS1.select(pol).toBands().reduceRegion('mean',ee_point,20)
        # Get the forest mask
        forestMask=get_forest_mask(date2,yearly_mask=yearly_mask, optical_mask=optical_mask, sar_mask=sar_mask, sar_tmp_mask=sar_tmp_mask)
        forestMaskBuffer=forestMask.focal_min(forestBuffer,"circle","meters")
        ###########################
        ### Stabilize S1 series ###
        ###########################
        # Stabilization 1: spatial p95 (Reiche, 2017)
        # def computep95(img):
        #     p95=img.reduceNeighborhood(
        #         reducer=ee.Reducer.percentile([95]),\
        #         kernel=ee.Kernel.square(searchDistance,'meters',False),\
        #         optimization='window')
        #     return p95.copyProperties(img,['system:time_start'])
        # p95col=colS1.select(pol).map(computep95)
        # #p95mean=p95col.mean()
        # #p95depart=p95col.map(lambda img:img.subtract(p95mean).copyProperties(img,['system:time_start']))
        # colS1_stable1=colS1.combine(p95col).map(lambda img: img.select(0).divide(img.select(1)).copyProperties(img,['system:time_start']))
        # Stabilization 2: 1st order harmonic fitting (Reiche, 2018)
        def addVariables(image):
            years = image.date().difference(date1, 'year')
            return image.addBands(ee.Image(years).rename('t')).addBands(ee.Image.constant(1)).float()
        colS1_st2 = colS1.select(pol).map(addVariables)
        #harmonicIndependents = ee.List(['constant', 't', 'cos', 'sin']) # Won't use trend
        harmonicIndependents = ee.List(['constant', 'cos', 'sin'])
        dependent = pol
        def addHarmonizedVars(image):
            timeRadians = image.select('t').multiply(2 * math.pi)
            return image.addBands(timeRadians.cos().rename('cos')).addBands(timeRadians.sin().rename('sin'))
        harmonicsS1=colS1_st2.map(addHarmonizedVars)
        harmonicLR = harmonicsS1.filterDate(date1,date2)\
            .select(harmonicIndependents.add(dependent))\
            .reduce(ee.Reducer.linearRegression(harmonicIndependents.length(), 1))
        # Turn the array image into a multi-band image of coefficients.
        harmonicLRCoefficients = harmonicLR.select('coefficients')\
            .arrayProject([0])\
            .arrayFlatten([harmonicIndependents])
        def computeFittedHarmonicNoIntercept(image):
            # fit = image.select(harmonicIndependents)\
            #         .multiply(harmonicLRCoefficients)\
            fit = image.select(['cos','sin'])\
                    .multiply(harmonicLRCoefficients.select(['cos','sin']))\
                    .reduce('sum')\
                    .rename('fitted')
            return image.addBands(fit)\
              .addBands(image.select(pol).subtract(fit).rename("residual"))
        fittedHarmonic=harmonicsS1.map(computeFittedHarmonicNoIntercept)
        colS1_stable2=fittedHarmonic.select("residual")
        # Stabilization 3: proposed here
        colS1_stable3=stabilize_SAR_time_series(colS1.select(pol),forestMaskBuffer,searchDistance)
        # COMPUTE POINT TIME SERIES
        #stabilized_time_series1=colS1_stable1.toBands().reduceRegion('mean',ee_point,20)
        stabilized_time_series2=colS1_stable2.toBands().reduceRegion('mean',ee_point,20)
        stabilized_time_series3=colS1_stable3.toBands().reduceRegion('mean',ee_point,20)
        ###### END STABILIZATION ##############
        # SPECKLE FILTERING
        median5=makeMedianFilter(5)
        colS1_f=QueganYuFilter(colS1.select(pol),median5).map(refinedLeeFilter)
        colS1_stable2_f=QueganYuFilter(colS1_stable2,median5).map(refinedLeeFilter)
        colS1_stable3_f=QueganYuFilter(colS1_stable3,median5).map(refinedLeeFilter)
        
        original_time_series_f=colS1_f.toBands().reduceRegion('mean',ee_point,20)
        stabilized_time_series2_f=colS1_stable2_f.toBands().reduceRegion('mean',ee_point,20)
        stabilized_time_series3_f=colS1_stable3_f.toBands().reduceRegion('mean',ee_point,20)
    else:
        original_time_series=\
        original_time_series_f=\
        stabilized_time_series2=\
        stabilized_time_series3=\
        stabilized_time_series2_f=\
        stabilized_time_series3_f=None    
    return LIA,original_time_series,original_time_series_f,\
                stabilized_time_series2,\
                stabilized_time_series3,\
                stabilized_time_series2_f,\
                stabilized_time_series3_f    
                
def get_point_timeseries_filter_2stab(pol,ee_point,ee_date_detect,yearly_mask,optical_mask,sar_mask,sar_tmp_mask,\
                        detection_threshold=3,generalScale=20,forestBuffer=100,searchDistance=5000,\
                        detectPeriodMonths=3,learningPeriodYears=3,stabilize=1,\
                        stabilize_scale=0,shrink=0,include_S1B=0,\
                        adaptative_threshold=None):
    logging.getLogger('googleapicliet.discovery_cache').setLevel(logging.ERROR)
    AOI=ee_point.buffer(10000)
    # Define relevant dates
    date2=ee_date_detect
    date3=date2.advance(detectPeriodMonths,'months')# Begin of detection (inclusive)
    date1=date2.advance(-1*learningPeriodYears,'years') # Begin of training colection (inclusive)
    date1str=date1.format('YYYY-MM-dd').getInfo()
    date2str=date2.format('YYYY-MM-dd').getInfo()
    date3str=date3.format('YYYY-MM-dd').getInfo()
    # First status messages
    print ("Processing point")
    print ("Learning period start: "+date1str)
    print ("Detection start: "+date2str)
    print ("Detection end: "+date3str)
    # Define S1 colection from EE
    colS1=getS1dataFloat(date1,date3) \
          .filterBounds(AOI)\
          .sort('system:time_start')\
          .map(lambda img: img.updateMask(img.gt(.00002))) \
          .map(toGamma0natural) \
          .select(["VHg0","VVg0","LIA"])
          
    LIA=colS1.select(["LIA"]).mean().reduceRegion('mean',ee_point,generalScale).getInfo().get("LIA")
    print(f"Local Incidence Angle: {LIA}")
    if not include_S1B:
        colS1=colS1.filterMetadata("platform_number","equals","A")
    # Leave only the most represented orbit in every pixel
    # add relative orbit number
    def addOrbitBand(img):
        orb=ee.Number(img.get('relativeOrbitNumber_start'))
        return img.addBands(ee.Image(orb).rename('Orbit').toInt())
    colS1=colS1.map(addOrbitBand)
    orbit_mode=colS1.select('Orbit').reduce(ee.Reducer.mode())
    colS1=colS1.map(lambda img: img.updateMask(img.select('Orbit').eq(orbit_mode))).select(["VHg0","VVg0","LIA"])
    colSize=colS1.size().getInfo()    
    print ("# of collected images: "+str(colSize))
    if colSize>0:
        original_time_series=colS1.select(pol).toBands().reduceRegion('mean',ee_point,20)
        # SPECKLE FILTERING
        median5=makeMedianFilter(5)
        def renameBand(img):
            return img.rename(pol)
        colS1_f=QueganYuFilter(colS1.select(pol),median5).map(refinedLeeFilter).map(renameBand)
        # Get the forest mask
        forestMask=get_forest_mask(date2,yearly_mask=yearly_mask, optical_mask=optical_mask, sar_mask=sar_mask, sar_tmp_mask=sar_tmp_mask)
        forestMaskBuffer=forestMask.focal_min(forestBuffer,"circle","meters")
        ###########################
        ### Stabilize S1 series ###
        ###########################
        def addVariables(image):
            years = image.date().difference(date1, 'year')
            return image.addBands(ee.Image(years).rename('t')).addBands(ee.Image.constant(1)).float()
        colS1_st2 = colS1_f.select(pol).map(addVariables)
        #harmonicIndependents = ee.List(['constant', 't', 'cos', 'sin']) # Won't use trend
        harmonicIndependents = ee.List(['constant', 'cos', 'sin'])
        dependent = pol
        def addHarmonizedVars(image):
            timeRadians = image.select('t').multiply(2 * math.pi)
            return image.addBands(timeRadians.cos().rename('cos')).addBands(timeRadians.sin().rename('sin'))
        harmonicsS1=colS1_st2.map(addHarmonizedVars)
        harmonicLR = harmonicsS1.filterDate(date1,date2)\
            .select(harmonicIndependents.add(dependent))\
            .reduce(ee.Reducer.linearRegression(harmonicIndependents.length(), 1))
        # Turn the array image into a multi-band image of coefficients.
        harmonicLRCoefficients = harmonicLR.select('coefficients')\
            .arrayProject([0])\
            .arrayFlatten([harmonicIndependents])
        def computeFittedHarmonicNoIntercept(image):
            # fit = image.select(harmonicIndependents)\
            #         .multiply(harmonicLRCoefficients)\
            fit = image.select(['cos','sin'])\
                    .multiply(harmonicLRCoefficients.select(['cos','sin']))\
                    .reduce('sum')\
                    .rename('fitted')
            return image.addBands(fit)\
              .addBands(image.select(pol).subtract(fit).rename("residual"))
        fittedHarmonic=harmonicsS1.map(computeFittedHarmonicNoIntercept)
        colS1_stable2=fittedHarmonic.select("residual")
        # Stabilization 3: proposed here
        colS1_stable3=stabilize_SAR_time_series(colS1_f.select(pol),forestMaskBuffer,searchDistance)
        # COMPUTE POINT TIME SERIES
        #stabilized_time_series1=colS1_stable1.toBands().reduceRegion('mean',ee_point,20)
        stabilized_time_series2=colS1_stable2.toBands().reduceRegion('mean',ee_point,20)
        stabilized_time_series3=colS1_stable3.toBands().reduceRegion('mean',ee_point,20)
        ###### END STABILIZATION ##############

        
        original_time_series_f=colS1_f.toBands().reduceRegion('mean',ee_point,20)
        stabilized_time_series2_f=colS1_stable2.toBands().reduceRegion('mean',ee_point,20)
        stabilized_time_series3_f=colS1_stable3.toBands().reduceRegion('mean',ee_point,20)
    else:
        original_time_series=\
        original_time_series_f=\
        stabilized_time_series2_f=\
        stabilized_time_series3_f=None    
    return LIA,original_time_series,original_time_series_f,\
                stabilized_time_series2_f,\
                stabilized_time_series3_f    
                
                

    # FILTERING
    #median5=makeMedianFilter(5)
    #colS1_f=QueganYuFilter(colS1_stable2,median5).map(refinedLeeFilter)
    #filtered_time_series=colS1_f.toBands().reduceRegion('mean',ee_point,20)
    # add date band and filters out extreme Local Incidence Angles to avoid misdetections
    #colS1_f2=colS1_f.combine((colS1).select("LIA")) \
    #    .map(lambda img: img.updateMask(img.select("LIA").gt(28)).updateMask(img.select("LIA").lt(50))).select(0)\
    #    .map(lambda img: img.updateMask(forestMask).addBands(ee.Image(img.date().difference("2020-01-01","days")).int16().rename("julian_day")))
#    return LIA,original_time_series,stabilized_time_series,filtered_time_series,loglikelihood_time_series
    

    # # DETECTION
    # # fix Forest and Non Forest distribuion parameters using temporal filtered collection
    # print ("Initializing detection")
    # learnCol=colS1_f2.select(0).filterDate(date1,date2)
    # nparamsF=learnCol.reduce(ee.Reducer.mean().combine(reducer2=ee.Reducer.stdDev(),sharedInputs= True))
    # nparamsNF=nparamsF.select(0).multiply(0.787952).subtract(0.013476).addBands(.005378298)
    # # computes loglikehood for every image of the detection period using spatio-temporal collection
    # def get_likelihood(img):
    #   pF=getNormalDistPdf(img,nparamsF)
    #   pNF=getNormalDistPdf(img,nparamsNF)
    #   pF=pF.where(img.gt(nparamsF.select(0)),ee.Image(0.39894228).divide(nparamsF.select(1)))
    #   pNF=pNF.where(img.lt(nparamsNF.select(0)),ee.Image(0.39894228).divide(nparamsNF.select(1)))
    #   NF_loglikelihood=pNF.divide(pF).log()
    #   return NF_loglikelihood.rename("NF_loglikelihood").copyProperties(img,['system:time_start'])
    # colS1_f_NFloglikelihood=colS1_f.select(0).filterDate(date1,date3).map(get_likelihood)
    # loglikelihood_time_series=colS1_f_NFloglikelihood.toBands().reduceRegion('mean',ee_point,20)
    
    # # creates warning raster  
    # if adaptative_threshold==None:
    #     detection_threshold=ee.Image(detection_threshold)
    # else:
    #     print ("Using adaptative thresholding")
    #     detection_threshold=ee.Image(adaptative_threshold)
    # colS1_f_NFloglikelihood_boolean=colS1_f_NFloglikelihood.map(lambda img:img.gte(detection_threshold))
    # warningCount=colS1_f_NFloglikelihood_boolean.sum()
    # nimg=colS1_f_NFloglikelihood_boolean.toArray().arrayProject([0]).arrayArgmax().arrayFlatten([["n"]])
    # dd=ee.Image(colS1_f2.filterDate(date2,date3).select(['julian_day']).toArray().arrayGet(nimg.addBands(0)))
    # contractPixels=1
    # dilatePixels=1
    # opening=2
    # warning_raster=ee.Image(0).where(warningCount.gte(1),1)\
    #         .focal_min(contractPixels,"square").focal_max(dilatePixels,"square")\
    #         .focal_max(opening,"circle").focal_min(opening,"circle")\
    #         .updateMask(forestMask)\
    #         .rename("alert")\
    #   .addBands(warningCount.updateMask(warningCount).rename("n_alerts"))\
    #   .addBands(dd.updateMask(warningCount).rename("daydetec"))\
    #   .addBands(colS1_f_NFloglikelihood.max().updateMask(warningCount).rename("prob_max").multiply(10))
     
    # warning_raster=warning_raster.updateMask(warning_raster.select(0)).toInt16()
