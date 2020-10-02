var BLA = ee.FeatureCollection("users/juandb/PRODES2019/brazilian_legal_amazon");
var warnings_sampling_image = ee.Image('users/juandb/DETER_SAR_EXP/warnings_sampling_image_v2')
var generalScale=20
var seed=23532542
var n=3000
var samples=warnings_sampling_image.stratifiedSample({
  numPoints:n,
  classBand:'constant',
  region:BLA,
  seed:seed,
  scale:generalScale,
  geometries:true
  })
//print ("Samples",samples)
Map.addLayer(samples,{},'Samples')
Export.table.toAsset(samples, 'samples_deforested_n'+n, 'DETER_SAR_EXP/samples_deforested_n'+n)