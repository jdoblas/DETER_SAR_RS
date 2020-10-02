var forest_2018_inpe = ee.Image("users/juandb/PESQUISA_DETER/floresta_prodes_2018_compress"),
    static_no_forest_mask = ee.Image("users/juandb/DETER_SAR_ANCILLARY/no_forest_mask_2019_raster_edit_23jun20"),
    deter_may2020 = ee.Image("users/juandb/DETER_SAR_ANCILLARY/deter_amz_public_2020Aug27_CR_raster"),
    BLA = ee.FeatureCollection("users/juandb/PRODES2019/brazilian_legal_amazon");
var generalScale=20
var forest_median_ndvi_threshold=0.85
Map.addLayer(forest_2018_inpe)
Map.addLayer(static_no_forest_mask)
Map.addLayer(deter_may2020)
var emersed_forest=forest_2018_inpe
                    .updateMask(static_no_forest_mask.unmask(0).not())
                    .updateMask(deter_may2020.gt(0).unmask(0).not())
Map.addLayer(emersed_forest,{},"emersed forest")
var kernel_forests=emersed_forest.unmask(0).focal_min(5000, 'circle', 'meters').selfMask()
Map.addLayer(kernel_forests,{},"Kernel forest")

// Selecting only green pixels
var forest_median_ndvi_threshold=0.85
var l8_sr1 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")
var l8_sr2 = ee.ImageCollection("LANDSAT/LC08/C01/T2_SR")
var l8_sr=l8_sr1.merge(l8_sr2).filterDate('2019-01-01','2020-01-01')
var l8_sr_masked=l8_sr.map(function(img){return img.updateMask(img.select(['pixel_qa']).eq(322))})
var l8_sr_masked_ndvi_median=l8_sr_masked.map(function(img){return img.normalizedDifference(['B5','B4'])
  .copyProperties(img,['system:time_start'])}).median()
Map.addLayer(l8_sr_masked_ndvi_median,{min:0.4,max:1,palette:['red','white','green']},'l8_sr_masked_ndvi_median')

var kernel_forests_green=kernel_forests.addBands(l8_sr_masked_ndvi_median.rename('median_NDVI')).updateMask(l8_sr_masked_ndvi_median.gt(forest_median_ndvi_threshold))
var scale=1000
var a=kernel_forests_green.select(0).reduceRegion(ee.Reducer.sum(),BLA,scale)
print(a)
print ("total kernel forest area",ee.Number(a.values().get(0)).multiply(scale^2))
Map.addLayer(kernel_forests_green,{bands:['b1'],palette:['green']},"kernel_forests_green")
Export.image.toAsset({
  image: kernel_forests_green,
  region: BLA,
  scale: generalScale,
  description: 'invariant_forest_sampling_image',
  assetId: 'DETER_SAR_EXP/invariant_forest_sampling_image',
  maxPixels: 216967828660
})
