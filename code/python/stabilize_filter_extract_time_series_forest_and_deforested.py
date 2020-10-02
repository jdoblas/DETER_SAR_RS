# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 11:44:48 2020

@author: juand
"""
import ee,sys,datetime
from get_time_series_stab_harmonic_spatialmean_filter import get_point_timeseries_2stab_filter

USAGE = f"Usage: python {sys.argv[0]} point_asset_id is_forest learning_period band output_prefix init_point end_point"

args = sys.argv[1:]
if not args:
    raise SystemExit(USAGE)
if len(args)<5:
    print(USAGE)

stabilize=1
include_S1B=0
points_collection=args[0]
is_forest=int(args[1])
learningPeriodYears=int(args[2])
pol=args[3]
output_prefix=args[4]

yearly_mask_asset="users/juandb/DETER_SAR_ANCILLARY/no_forest_mask_2019_raster_edit_21jun20"
optical_mask_asset="users/juandb/DETER_SAR_ANCILLARY/deter_amz_public_2020Aug27_CR_raster"
sar_mask_asset="users/juandb/DETER_SAR_ANCILLARY/DETERSAR_CR2_MASK_raster"
# test point 
ee.Initialize()
time_start=datetime.datetime.now()
npontos=ee.FeatureCollection(points_collection).size().getInfo()
if len(args)>5:
    init_point=args[5]
else:
    init_point=0
if len(args)==7:
    end_point=args[6]
else:
    end_point=npontos
if is_forest: 
    lista_pontos=ee.FeatureCollection(points_collection).randomColumn().toList(npontos)
else:
    lista_pontos = ee.FeatureCollection(points_collection).toList(npontos)
    
if is_forest:
    output_name=output_prefix+"_forestinvariant_"+str(learningPeriodYears)+"_"+pol+".csv"
else:
    output_name=output_prefix+"_deforested_"+str(learningPeriodYears)+"_"+pol+".csv"
#### EXTRACTION BEGINS HERE ####    
print (f"Computing stabilized {pol} TS with {learningPeriodYears}yr period")
with open(output_name,"a") as f:
    for j in range(init_point,end_point):
        point = ee.Feature(lista_pontos.get(j))
        print (f"Point {j}")
        if is_forest:
            ee_begin_detect_date=ee.Date('2019-01-01').advance(ee.Number(point.get('random')).multiply(365),'days') # We generate a random deforestation date
        else:
            ee_begin_detect_date=ee.Date(point.get('before_dt'))
        LIA,original,original_f,stabilized2,stabilized3,filtered2,filtered3=get_point_timeseries_2stab_filter(pol,point.geometry(),ee_begin_detect_date,\
                                                learningPeriodYears=learningPeriodYears,\
                                                detectPeriodMonths=12,\
                                                stabilize=stabilize,\
                                                yearly_mask=yearly_mask_asset,\
                                                optical_mask=optical_mask_asset,\
                                                sar_mask=sar_mask_asset,\
                                                sar_tmp_mask=ee.FeatureCollection([]))
        if original!=None:
            a=original.getInfo()
            af=original_f.getInfo()
            b=stabilized2.getInfo()
            c=stabilized3.getInfo()
            d=filtered2.getInfo()
            e=filtered3.getInfo()
            values=list(a.keys())
            anp=list(a.values())
            afnp=list(af.values())
            bnp=list(b.values())
            cnp=list(c.values())
            dnp=list(d.values())
            enp=list(e.values())
            if is_forest:
                write_vector=[j,LIA,point.geometry().coordinates().getInfo(),\
                                point.get('median_NDVI').getInfo(),\
                                ee_begin_detect_date.millis().getInfo(),\
                                values,anp,afnp,bnp,cnp,dnp,enp]
            else:
                write_vector=[j,LIA,point.geometry().coordinates().getInfo(),\
                                point.get('NDVI_after').getInfo(),\
                                point.get('NDVI_before').getInfo(),\
                                point.get('after_dt').getInfo(),\
                                point.get('before_dt').getInfo(),\
                                point.get('detection_dt').getInfo(),\
                                values,anp,afnp,bnp,cnp,dnp,enp]
                
            f.write(str(write_vector)+'\n')
            f.flush()
time_end=datetime.datetime.now()
print ("elapsed time ",time_end-time_start)
