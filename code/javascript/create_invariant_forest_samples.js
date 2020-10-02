var BLA = ee.FeatureCollection("users/juandb/PRODES2019/brazilian_legal_amazon");
var warnings_sampling_image = ee.Image('users/detersaree/invariant_forest_sampling_image_kernel_1km')
var generalScale=20
var seed=12312
var n=3000
print (warnings_sampling_image)
Map.addLayer(warnings_sampling_image)
var samples=warnings_sampling_image.stratifiedSample({
  numPoints:n,
  classBand:'b1',
  region:BLA,
  seed:seed,
  scale:generalScale,
  geometries:true
  })
//print ("Samples",samples)
Map.addLayer(samples,{},'Samples')
Export.table.toAsset(samples, 'samples_invariant_n'+n, 'DETER_SAR_EXP/samples_invariant_n'+n)