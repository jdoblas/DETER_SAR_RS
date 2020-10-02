# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 09:38:21 2020

@author: juand
"""
import ee,datetime

def get_forest_mask(date_start_detection,yearly_mask,optical_mask,sar_mask,sar_tmp_mask):
    """
    Returns a forest mask, computed as the negation of the addition of all available rasterized non-forest masks:
        - Acummulated 2019 deforestation (PRODES/INPE), hidrology (INPE), non forest (INPE) and varzea mask (INPE)
        - Ongoing optical deforestation warnings (DETER), until date_start_detection
        - DETERSAR CR2 warnings until the actual run
        - DETERSAR_CR2 warnings computed on the actual run 
    Args:
        date_start_detection (ee.Date): Date of the start of SAR detection. Used to mask out any posterior optical deforestation data.
        yearly_mask
        optical_mask (STR): Current rasterized mask of DETER. Should have unique band representing the number of elapsed days from 2018-01-01.
                            Use 'precompute_deforest_mask' to compute this image from DETER polygons.
        sar_mask (STR):    Accumulated detections of DETERSAR
        sar_tmp_mask (ee.FeatureCollection): Feature Collection of the CR2 polygons detected during the current run
    Returns:
        ee.Image: mask representing forested pixels.

    """
    deforestationMaskRaster=make1sFromMask(ee.Image(yearly_mask)).unmask(0,False)
    days_from_1jan2018=date_start_detection.difference(ee.Date('2018-01-01'),'day')
    if optical_mask==None:
        deforestationMaskRaster_DETER=ee.Image(0)
    else:
        deforestationMaskRaster_DETER=ee.Image(optical_mask).lt(ee.Image(days_from_1jan2018)).unmask(0,False)
    if sar_mask==None:
        deforestationMaskRaster_DETERSAR=ee.Image(0)
    else:
        deforestationMaskRaster_DETERSAR=make1sFromMask(ee.Image(sar_mask)).unmask(0,False)
    deforestationMaskRaster_DETERSAR_TMP=ee.Image(sar_tmp_mask.map(lambda ft:ft.set('desm',1)).reduceToImage(['desm'],'first')).unmask(0,False)
    deforestationMask=deforestationMaskRaster.Or(deforestationMaskRaster_DETER).Or(deforestationMaskRaster_DETERSAR).Or(deforestationMaskRaster_DETERSAR_TMP)
    return deforestationMask.Not()

def make1sFromMask(img):
    return img.selfMask().add(1).divide(img.add(1))

def precompute_deforest_mask(warning_vector_asset="users/juandb/DETER_SAR_PRODUCAO/deter-amz-public-2020Mai26",output_asset="users/juandb/DETER_SAR_PRODUCAO/deter_amz_public_2020Mai26_CR_raster"):
    DETER_CR = ee.FeatureCollection(warning_vector_asset) \
                .filterMetadata("CLASSNAME","equals","DESMATAMENTO_CR")
    DETER_CR_VEG = ee.FeatureCollection(warning_vector_asset) \
                .filterMetadata("CLASSNAME","equals","DESMATAMENTO_VEG")
    def add_date(ft):return ft.set("days_from_1jan2018",ee.Date(ft.get("VIEW_DATE")).difference(ee.Date('2018-01-01'),'day'))                  
    DETER=DETER_CR.merge(DETER_CR_VEG).map(add_date).filterMetadata("days_from_1jan2018","greater_than",212)
    print (DETER.size().getInfo())
    deforestationMaskRaster_DETER=DETER.reduceToImage(["days_from_1jan2018"],ee.Reducer.mean()).rename('days_from_2018').toInt16()
    print(deforestationMaskRaster_DETER.getInfo())
    task=ee.batch.Export.image.toAsset(image=deforestationMaskRaster_DETER,description="export_forest_mask",assetId=output_asset,scale=20,maxPixels=3400123231231)
    return task

def rasterize_deforestation_polygons(pols_asset="users/juandb/DETER_SAR_PRODUCAO/DETERSAR_CR2_TRANSMITIDOS_4AREAS_ACUMULADO",output_asset="users/juandb/DETER_SAR_PRODUCAO/DETERSAR_CR2_TRANSMITIDOS_4AREAS_ACUMULADO_raster"):
    raster=ee.FeatureCollection(pols_asset).map(lambda ft:ft.set('desm',1)).reduceToImage(['desm'],'first')
    task=ee.batch.Export.image.toAsset(image=raster,description="rasterize",assetId=output_asset,scale=20,maxPixels=3400123231231)
    return task

if __name__== '__main__':
    task16=rasterize_deforestation_polygons()
    