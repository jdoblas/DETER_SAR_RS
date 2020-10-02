# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 10:27:50 2020

@author: juand
"""

import ee,math

def stabilize_SAR_time_series (collection, mask, search_distance_meters):
    # collection should have only one band of linear scale backscattering
    def compute_forest_mean(img):
        forest_mean=img.updateMask(mask).reduceNeighborhood(
          reducer = ee.Reducer.mean(),
          kernel = ee.Kernel.square(search_distance_meters,"meters",False),
          inputWeight = "mask",
          skipMasked = False,
          optimization = "boxcar"
          ).copyProperties(img,["system:time_start"])
        return img.addBands(ee.Image(forest_mean))
    colS1Forest=collection.map(compute_forest_mean)
    def compute_coef(img):
      return img.addBands(img.select(1).divide(colS1Forest.select(1).mean()).unmask(1).copyProperties(img,["system:time_start"]))
    colS1Forest=colS1Forest.map(compute_coef)
       
    def compute_corrected_backscattering(img):
      return img.select(0).divide(img.select(2)).copyProperties(img,["system:time_start"])
    return colS1Forest.map(compute_corrected_backscattering)

def stabilize_SAR_time_series_upscale (collection, mask, search_distance_meters,scale):
    # collection should have only one band of linear scale backscattering
    def compute_forest_mean(img):
          forest_mean=img.updateMask(mask).reduceNeighborhood(
          reducer = ee.Reducer.mean(),
          kernel = ee.Kernel.square(search_distance_meters,"meters",False),
          inputWeight = "mask",
          skipMasked = False,
          optimization = "boxcar"
          ).reproject(crs=img.projection(),scale=scale).copyProperties(img,["system:time_start"])
          #).reduceResolution('mean', True, 64).reproject(crs=img.projection(),scale=scale).copyProperties(img,["system:time_start"])
          return img.addBands(ee.Image(forest_mean))
    colS1Forest=collection.map(compute_forest_mean)
    def compute_coef(img):
      return img.addBands(img.select(1).divide(colS1Forest.select(1).mean()).unmask(1).copyProperties(img,["system:time_start"]))
    colS1Forest=colS1Forest.map(compute_coef)
       
    def compute_corrected_backscattering(img):
      return img.select(0).divide(img.select(2)).copyProperties(img,["system:time_start"])
    return colS1Forest.map(compute_corrected_backscattering)

def stabilize_SAR_time_series_harmonic (collection,date1,date2,pol):
    def addVariables(image):
        years = image.date().difference(date1, 'year')
        return image.addBands(ee.Image(years).rename('t')).addBands(ee.Image.constant(1)).float()
    colS1_st2 = collection.map(addVariables)
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
    return fittedHarmonic.select("residual")
    