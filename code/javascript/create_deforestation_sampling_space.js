var BLA = ee.FeatureCollection("users/juandb/PRODES2019/brazilian_legal_amazon");
var generalScale=20
var warnings=ee.FeatureCollection('users/juandb/DETER_SAR_EXP/mapbiomas_alerts_allattribs_2019_AMAZONIA')
//print (warnings.first())
// Selects 2019 alerts, alerts > 1ha, have less than 4 months image interval, 
// date of forested images is before date of detection
warnings=warnings.map(function(ft){
                        ft=ft.set('detection_dt',ee.Date(ft.get('DataDetec')).millis())
                        var date_before=ft.get('before_dt')
                        var date_detect=ft.get('detection_dt')
                        var dif=ee.Number(date_detect).subtract(ee.Number(date_before))
                        return ft.set('dif_detect',dif)
                      })
print (warnings.first())
var filtered_warnings=warnings
                      .filterMetadata('AnoDetec','equals',2019)
                      .filterMetadata('AreaHa','greater_than',1)
                      .filterMetadata('days_inter','less_than',121)
                      .filterMetadata('after_dt','not_equals',null)
                      .filterMetadata('before_dt','not_equals',null)
                      .filterMetadata('dif_detect','greater_than',0)
Map.addLayer(filtered_warnings,{},"Filtered warning polygons")
print (filtered_warnings.size())
// Samples 10000 points
var seed=234234
//var initial_samples=ee.FeatureCollection.randomPoints(filtered_warnings.geometry(),10000,seed)
//Map.addLayer(initial_samples)
// Filters out points that weren't deforested after the NDVI
var ndvi_forest_threshold=.82
var ndvi_nonforest_threshold=.69
var l8_sr1 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")
var l8_sr2 = ee.ImageCollection("LANDSAT/LC08/C01/T2_SR")

var ft=filtered_warnings.first()
print (ft)
var mask_collection=ee.ImageCollection(filtered_warnings.map(function(ft){
  //print (ft)
  //Map.centerObject(ft)
  var date1=ee.Date(ft.get('before_dt'))
  var date2=ee.Date(ft.get('after_dt'))
  //print (date1)
  //print (date2)
  var l8_sr=l8_sr1.merge(l8_sr2)
  var l8_sr_masked=l8_sr.map(function(img){return img.updateMask(img.select(['pixel_qa']).eq(322))})
  var l8_sr_masked_ndvi=l8_sr_masked.map(function(img){return img.normalizedDifference(['B5','B4'])
    .copyProperties(img,['system:time_start'])
  })
  //print (l8_sr_masked_ndvi.first())
  var ndvi_before=l8_sr_masked_ndvi.filterDate(date1.advance(-3,'months'),date1).median()
  var ndvi_after=l8_sr_masked_ndvi.filterDate(date2,date2.advance(3,'months')).median()
  //Map.addLayer(ndvi_before,{min:0.4,max:1,palette:['red','white','green']},'NDVI before')
  //Map.addLayer(ndvi_after,{min:0.4,max:1,palette:['red','white','green']},'NDVI after')
  var mask=ndvi_before.gte(ndvi_forest_threshold).multiply(ndvi_after.lt(ndvi_nonforest_threshold))
              .addBands(ndvi_before.rename("NDVI_before"))
              .addBands(ndvi_after.rename("NDVI_after"))
              .clip(ft)
  return mask
})).mosaic()
var mask=ee.Image(0).blend(mask_collection.select(0)).focal_min(20,'circle', 'meters').selfMask()
Map.addLayer(mask,{},"Mascara NDVI")
// Rasterize filtered_warnings
var properties=['detection_dt','before_dt', 'after_dt']
var reducer=ee.Reducer.first().forEach(properties)
var filtered_warnings_raster=filtered_warnings
                              .reduceToImage(properties,reducer)
                              .addBands(mask_collection.select(1,2))
                              .addBands(ee.Image(1))
                              .updateMask(mask)
                              
Map.addLayer(filtered_warnings_raster,{},"Rasterized warnings")
Export.image.toAsset({
  image: filtered_warnings_raster,
  region: BLA,
  scale: generalScale,
  description: 'warnings_sampling_image',
  assetId: 'DETER_SAR_EXP/warnings_sampling_image',
  maxPixels: 216967828660
})
