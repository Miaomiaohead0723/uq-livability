var roi = table
var year='2015'
var st_date=year+'-01-01'
var en_date=year+'-12-31'
var scale_my=500
var kernelSize = 10
var filter=ee.Filter.and(
                ee.Filter.date(st_date,en_date),
                ee.Filter.bounds(roi)
             )
Map.centerObject(roi, 6);

var visRate = {
  min: 0,
  max: 1,
  palette: [
    '000080', '0000d9', '4000ff', '8000ff', '0080ff', '00ffff', '00ff80',
    '80ff00', 'daff00', 'ffff00', 'fff500', 'ffda00', 'ffb000', 'ffa400',
    'ff4f00', 'ff2500', 'ff0a00', 'ff00ff'
  ]
};
// print(Number(Math.random()))

function MonteCarlo(idx){
  var idx = ee.Number(idx)
  // print(idx)
  var ran1 = ee.Number(randnum1.get(idx))
  var ran2 = ee.Number(randnum2.get(idx))
  var ran3 = ee.Number(randnum3.get(idx))
  var ran4 = ee.Number(randnum4.get(idx))
  var ran5 = ee.Number(randnum5.get(idx))
  var ran6 = ee.Number(randnum6.get(idx))
  var ran7 = ee.Number(randnum7.get(idx))
  var ran8 = ee.Number(randnum8.get(idx))
  var ran9 = ee.Number(randnum9.get(idx))

  // print(ran1,ran2,ran3)

/************************************************************************
 * 
 *                           dataset 
 *  
 * *********************************************************************/
var SRTM = ee.Image('USGS/SRTMGL1_003').multiply(ran1)
var ERA5 = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR").filter(filter).median().clip(roi).multiply(ran2);
var sand = ee.Image("projects/soilgrids-isric/sand_mean").clip(roi).select("sand_0-5cm_mean").divide(10).unmask(0).multiply(ran3);
var silt = ee.Image("projects/soilgrids-isric/silt_mean").clip(roi).select("silt_0-5cm_mean").divide(10).unmask(0).multiply(ran4);
var clay = ee.Image("projects/soilgrids-isric/clay_mean").clip(roi).select("clay_0-5cm_mean").divide(10).unmask(0).multiply(ran5);
var col_a1 = ee.ImageCollection('MODIS/006/MOD09A1').filter(filter);  // wet ndvi ndbsi
var al_single=col_a1.median().select(['sur_refl_b01','sur_refl_b02','sur_refl_b03','sur_refl_b04','sur_refl_b05','sur_refl_b06','sur_refl_b07'],
                                     ["red","nir1","blue","green","nir2","swir1","swir2"]).clip(roi).multiply(ran6);
var pop = ee.ImageCollection('WorldPop/GP/100m/pop').filter(ee.Filter.date('2015-01-01', '2015-12-31')).median().clip(roi).multiply(ran7);
var gdp_ppp = ee.Image("projects/sat-io/open-datasets/GRIDDED_HDI_GDP/GDP_PPP_1990_2015_5arcmin_v2").multiply(ran8);
var era5_tp = ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD')
                  .select('precipitation')
                  .filter(filter)
                  .sum().multiply(10)
                  .clip(roi).multiply(ran9);
/************************************************************************
 * 
 *                           ecological environmental quality
 *  
 * *********************************************************************/

var a1_index={
   // mNDWI
   mNDWI:function(img){
      var mndwi = img.normalizedDifference(["green","swir1"]);  
      return img.addBands(mndwi.rename("mNDWI"));
   },
  // NDWI
   NDWI:function(img){
      var mndwi = img.normalizedDifference(["green","nir1"]);  
      return img.addBands(mndwi.rename("NDWI"));
   },
   // NDVI
   NDVI:function(img){
      var ndvi = img.normalizedDifference(["nir1","red"]);  //0.85 - 0.88 µm
      return img.addBands(ndvi.rename("NDVI"));
   },
   // WET
   WET:function WET(img){
        var wet =img.expression(
          "float(0.1147*b1+ 0.2489*b2+ 0.2408*b3+ 0.3132*b4- 0.3122*b5- 0.6416*b6- 0.5087*b7)",
        {
          "b1":img.select("red"),
          "b2":img.select("nir1"),
          "b3":img.select("blue"),
          "b4":img.select("green"),
          "b5":img.select("nir2"),
          "b6":img.select("swir1"),
          "b7":img.select("swir2"),
        }
        );
        return img.addBands(wet.multiply(0.0001).rename("WET"));
    },
    // NDBSI
    NDBSI:function NDBSI(img){ // NDBSI
      var ibi = img.expression(
        '((2*swir/(swir+nir))-((nir/(nir+red))+(green/(green+swir))))/((2*swir/(swir+nir))+((nir/(nir+red))+(green/(green+swir))))',
        {
          'swir':img.select('swir1'),//1.57 - 1.65 µm
          'nir':img.select('nir1'),
          'red':img.select('red'),
          'green':img.select('green')
        });
  
      var ndsi = img.expression(
        '((swir+red)-(blue+nir))/(swir+red+blue+nir)',
        {
          'swir':img.select('swir1'),
          'nir':img.select('nir1'),
          'red':img.select('red'),
          'green':img.select('green'),
          'blue':img.select('blue')
        });
        
      var ndbsi=ee.Image(ibi).add(ee.Image(ndsi)).divide(2)
      return img.addBands(ndsi.rename("NDBSI"));
    },
    
    MSAVI:function MSAVI(img){
      var msavi = img.expression(
        "(2*nir-1)",
        // "(2*nir-1)-pow(2*nir+1-8*(nir-red),0.5)",
        {
          'nir':img.select('nir1'),
          'red':img.select('red')
        });
      return img.addBands(msavi.rename("MSAVI"));
    },
    
    SI:function SI(img){
      var si = img.expression(
        "pow(blue*red,0.5)",
        {
          'blue':img.select('blue'),
          'red':img.select('red')
        });
      return img.addBands(si.rename("SI"));
    }
}


// Map.addLayer(al_single,{},'al_single')
var mndwi=a1_index.mNDWI(al_single).select("mNDWI");
var water_mask = mndwi.lt(0.2);

al_single = al_single.updateMask(water_mask);


var ndwi=a1_index.NDWI(al_single).select("NDWI")
var ndvi=a1_index.NDVI(al_single).select("NDVI")

var msavi=a1_index.MSAVI(al_single).select("MSAVI")
// Map.addLayer(msavi,{},'msavi')
var si=a1_index.SI(al_single).select("SI")

var bandNames = ['msavi', 'si'];
var combinedImage = msavi.addBands(si)
var soil = combinedImage.select([0, 1], bandNames);
// print('soil',soil)
// Map.addLayer(soil,{},'soil')
//标准化 
function normalization(image,region,scale){
   
   var num = image.reduceRegion({
        reducer:ee.Reducer.percentile([1,99]),
        geometry:region,
        scale:scale,
        maxPixels:1e13,
        // tileScale: 16
      })
  
// use unit scale to normalize the pixel values
  var unitScale = ee.ImageCollection.fromImages(
    image.bandNames().map(function(name){
    name = ee.String(name);
    var num_1 = ee.Number(num.get(name.cat('_p1'))).add(0.0001); // get the minimum cutoff value
    var num_99 = ee.Number(num.get(name.cat('_p99'))).add(0.0002); // get the maximum cutoff value
    var band = image.select(name);
    var max = num_99;
    var min = num_1;
    var band1=ee.Image(min).multiply(band.lt(min)).add(ee.Image(max).multiply(band.gt(max)))
                        .add(band.multiply(ee.Image(1).subtract(band.lt(min)).subtract(band.gt(max))))
    var result_band=band1.subtract(min).divide(max.subtract(min));
    return result_band;
  })).toBands().rename(image.bandNames()); // conver to band
    return unitScale;
}

var normalSoil=ee.Image(normalization(soil,roi,scale_my))

var expression = 'pow(pow(msavi-1,2)+pow(si,2),0.5)';

// Compute the expression for each pixel in the image or image collection
var msdi = normalSoil.expression(expression, {
  'msavi': normalSoil.select('msavi'),
  'si': normalSoil.select('si')
});
// print('normalSoil',normalSoil)
// print('msdi',msdi)
/************************************************************************
 * 
 *                           topography factor
 *  
 * *********************************************************************/

var elevation = SRTM.select('elevation').clip(roi).updateMask(water_mask).float();
var slope = ee.Terrain.slope(elevation);  //Slope
var aspect = ee.Terrain.aspect(elevation);  //Aspact

var ALT = elevation.reduceNeighborhood({
  reducer: ee.Reducer.mean(),
  kernel: ee.Kernel.square(kernelSize),
});

var maxH = elevation.reduceNeighborhood({
  reducer: ee.Reducer.max(),
  kernel: ee.Kernel.square(kernelSize),
});

var minH = elevation.reduceNeighborhood({
  reducer: ee.Reducer.min(),
  kernel: ee.Kernel.square(kernelSize),
});

var flat = slope.where(slope.lte(5), 1).where(slope.gt(5), 0);
var flatAreaRatio = flat.reduceNeighborhood({
  reducer: ee.Reducer.sum(),
  kernel: ee.Kernel.square(kernelSize),
});

var bandNames = ['ALT', 'maxH', 'minH', 'flatAreaRatio'];
var combinedImage = ALT.addBands(maxH).addBands(minH).addBands(flatAreaRatio)
var topography = combinedImage.select([0, 1, 2, 3], bandNames);

var expression = 'ALT / 100 + (maxH-minH)*(1-flatAreaRatio)/500';

var RDLS = topography.expression(expression, {
  'ALT': topography.select('ALT'),
  'maxH': topography.select('maxH'),
  'minH': topography.select('minH'),
  'flatAreaRatio': topography.select('flatAreaRatio')
});


/************************************************************************
 * 
 *                           climate factor
 *  
 * *********************************************************************/
var era5_2mt = ERA5.select('temperature_2m').subtract(273.15).updateMask(water_mask);
var era5_u_wind_10m = ERA5.select('u_component_of_wind_10m').updateMask(water_mask);
var era5_v_wind_10m = ERA5.select('v_component_of_wind_10m').updateMask(water_mask);

var expressionWCI = '13.12+0.6215*T-11.37*pow(pow(u,2)+pow(v,2), 0.08)+0.3965*T*pow(pow(u,2)+pow(v,2), 0.08)';

var WCI = ERA5.expression(expressionWCI, {
  'T': ERA5.select('temperature_2m').subtract(273.15),  // Replace 'ALT' with the actual band name for ALT
  'u': ERA5.select('u_component_of_wind_10m'),  // Replace 'maxH' with the actual band name for maxH
  'v': ERA5.select('v_component_of_wind_10m'),  // Replace 'minH' with the actual band name for minH
});

var expressionICHB = '1.8 * T + 32 - 0.55 * 0.5 * (1.8 * T - 26) - 3.2 * pow(pow(u,2)+pow(v,2), 0.5)';
var ICHB = ERA5.expression(expressionICHB, {
  'T': ERA5.select('temperature_2m').subtract(273.15),  // Replace 'ALT' with the actual band name for ALT
  'u': ERA5.select('u_component_of_wind_10m'),  // Replace 'maxH' with the actual band name for maxH
  'v': ERA5.select('v_component_of_wind_10m'),  // Replace 'minH' with the actual band name for minH
});


/************************************************************************
 * 
 *                           GDP
 *  
 * *********************************************************************/
pop = pop.updateMask(water_mask);
var GDP = gdp_ppp.select('b26').clip(roi).updateMask(water_mask);


/************************************************************************
 * 
 *                           dataset 
 *  
 * *********************************************************************/

                  var bandNames = ['elevation', 'slope'];
var combinedImage = elevation.addBands(slope)
var soil = combinedImage.select([0, 1], bandNames);



var slopeRank=slope.where(slope.lte(8),1)
                .where(slope.gt(8).and(slope.lte(15)),2)
                .where(slope.gt(15).and(slope.lte(25)),3)
                .where(slope.gt(25).and(slope.lte(35)),4)
                .where(slope.gt(35),5)
                
var elevationRank=elevation.where(elevation.lte(1000),1)
                .where(elevation.gt(1000).and(elevation.lte(1300)),2)
                .where(elevation.gt(1300).and(elevation.lte(1700)),3)
                .where(elevation.gt(1700).and(elevation.lte(2500)),4)
                .where(elevation.gt(2500),5)
                
var expressionSlopeLength = 'elevation / (sin(slope/180))';

var slopeLength = soil.expression(expressionSlopeLength, {
  'elevation': soil.select('elevation').multiply(Math.PI),
  'slope': soil.select('slope'), 
});

var slopeLengthRank = slopeLength.where(slopeLength.lte(150000),1)
                                  .where(slopeLength.gt(150000).and(slopeLength.lte(300000)),2)
                                  .where(slopeLength.gt(300000).and(slopeLength.lte(400000)),3)
                                  .where(slopeLength.gt(400000).and(slopeLength.lte(600000)),4)
                                  .where(slopeLength.gt(600000),5)


var blankData = ee.FeatureCollection([
  ee.Feature(roi, {})
]);
var filledData = blankData.map(function(feature) {
  return feature.set('value', 0);
});

var sandRank=sand.where(sand.lte(40),5)
                .where(sand.gt(40).and(sand.lte(50)),4)
                .where(sand.gt(50).and(sand.lte(60)),3)
                .where(sand.gt(60).and(sand.lte(80)),2)
                .where(sand.gt(80),1)
                
var siltRank=elevation.where(silt.lte(15),5)
                .where(silt.gt(15).and(silt),4)
                .where(silt.gt(25).and(silt.lte(30)),3)
                .where(silt.gt(30).and(silt.lte(35)),2)
                .where(silt.gt(35),1)
                
var clayRank=clay.where(clay.lte(10),5)
                .where(clay.gt(10).and(clay.lte(15)),4)
                .where(clay.gt(15).and(clay.lte(20)),3)
                .where(clay.gt(20).and(clay.lte(30)),2)
                .where(clay.gt(30),1)

var era5_tpRank=era5_tp.where(era5_tp.lte(2500),1)
                .where(era5_tp.gt(2500).and(era5_tp.lte(4000)),2)
                .where(era5_tp.gt(4000).and(era5_tp.lte(5000)),3)
                .where(era5_tp.gt(5000).and(era5_tp.lte(6500)),4)
                .where(era5_tp.gt(6500),1)


soil = soil.addBands(ndvi.rename('ndvi'))

var ndviRank=ndvi.where(ndvi.lte(0),5)
                .where(ndvi.gt(0).and(ndvi.lte(0.3)),4)
                .where(ndvi.gt(0.3).and(ndvi.lte(0.5)),3)
                .where(ndvi.gt(0.5),1)


var landUse = ee.ImageCollection('ESA/WorldCover/v100').first().clip(roi);

var landUseRank=landUse.where(landUse.eq(10),1)
                .where(landUse.eq(20),1)
                .where(landUse.eq(30),2)
                .where(landUse.eq(40),3)
                .where(landUse.eq(50),4)
                .where(landUse.eq(60),5)
                .where(landUse.eq(70),2)
                .where(landUse.eq(80),0)
                .where(landUse.eq(90),2)
                .where(landUse.eq(95),1)
                .where(landUse.eq(100),5)
                


var rankBandNames = ['slopeRank', 'elevationRank', 'slopeLengthRank','ndviRank','landUseRank','era5_tpRank','sandRank','siltRank','clayRank'];
var combinedImage = slopeRank.addBands(elevationRank)
                                  .addBands(slopeLengthRank)
                                  .addBands(ndviRank)
                                  .addBands(landUseRank)
                                  .addBands(era5_tpRank)
                                  .addBands(sandRank)
                                  .addBands(siltRank)
                                  .addBands(clayRank)
var rank = combinedImage.select([0,1,2,3,4,5,6,7,8], rankBandNames);

var expressionErosion = '0.15*slope+0.05*elevation+slopeLength*0.1+ndvi*0.15+landUse*0.1+tp*0.15+sand*0.1+silt*0.1+clay*0.1';

var erosion = rank.expression(expressionErosion, {
  'slope': rank.select('slopeRank'),
  'elevation': rank.select('elevationRank'),
  'slopeLength': rank.select('slopeLengthRank'),
  'ndvi': rank.select('ndviRank'),
  'landUse': rank.select('landUseRank'),
  'tp': rank.select('era5_tpRank'),
  'sand': rank.select('sandRank'),
  'silt': rank.select('siltRank'),
  'clay': rank.select('clayRank'),
});


/************************************************************************
 * 
 *                           show 
 *  
 * *********************************************************************/

//Map.addLayer(slopeRank, visParams, 'slopeRank');
//Map.addLayer(elevationRank, visParams, 'elevationRank');
//Map.addLayer(slopeLengthRank, visParams, 'slopeLengthRank');
//Map.addLayer(ndviRank, visParams, 'ndviRank');
//Map.addLayer(landUseRank, visParams, 'landUseRank');
//Map.addLayer(era5_tpRank, visParams, 'era5_tpRank');
//Map.addLayer(sandRank, visParams, 'sandRank');
//Map.addLayer(siltRank, visParams, 'siltRank');
//Map.addLayer(clayRank, visParams, 'clayRank');

//Map.addLayer(erosion, {min:0, max:5, palette: palette_rainbow}, 'erosion');

//Map.addLayer(slope, {min: 0, max: 60}, 'slope');
//Map.addLayer(aspect, {min: 0, max: 60}, 'aspect');
//Map.addLayer(ALT, visParam, 'ALT');
//Map.addLayer(maxH, visParam, 'maxH');
//Map.addLayer(minH, visParam, 'minH');
//Map.addLayer(flat, {min: 0, max: 1, palette: ['000000', 'FFFFFF']}, 'flat');
//Map.addLayer(flatAreaRatio, {min: 0, max: 1}, 'flatAreaRatio');
//Map.addLayer(RDLS, visParam2, 'RDLS');
//Map.addLayer(era5_2mt, vis2mt, 'era5_2mt');
//Map.addLayer(era5_u_wind_10m, visWind , 'era5_u_wind_10m');
//Map.addLayer(era5_v_wind_10m, visWind , 'era5_v_wind_10m');
//Map.addLayer(WCI, visWind,'WCI')
//Map.addLayer(ICHB, visICHB,'ICHB')
//Map.addLayer(mndwi, visRate,'mndwi')
//Map.addLayer(ndwi, visRate,'ndwi')
//Map.addLayer(ndvi, visRate,'ndvi')
//Map.addLayer(msavi, visParamMSDI,'msavi')
//Map.addLayer(si, visParamMSDI,'si')
//Map.addLayer(msdi, visParamMSDI,'msdi')
//Map.addLayer(pop, visPop, 'Population');
//Map.addLayer(GDP,{min:-347969610,max:411755463,palette: visGDP.extra.orngblue},'GDP')



//var resampledPop = pop.resample('bilinear').reproject({
//  crs: pop.projection(),
//  scale: 30
//});
//var resampledGDP = GDP.resample('bilinear').reproject({
//  crs: GDP.projection(),
//  scale: 30
//});


var normalErosion = ee.Image(normalization(erosion,roi,scale_my))
var normalRDLS = ee.Image(normalization(RDLS,roi,scale_my))
var normalWCI = ee.Image(normalization(WCI,roi,scale_my))
var normalICHB = ee.Image(normalization(ICHB,roi,scale_my)).subtract(1).multiply(-1)

var normalmsdi = ee.Image(normalization(msdi,roi,scale_my)).subtract(1).multiply(-1)
var normalndvi = ee.Image(normalization(ndvi,roi,scale_my))
var normalpop = ee.Image(normalization(pop,roi,scale_my))
var normalGDP = ee.Image(normalization(GDP,roi,scale_my))

var normalBandNames = ['normalErosion', 'normalRDLS', 'normalWCI','normalICHB','normalmsdi','normalndvi','normalpop','normalGDP'];
var combinedImage = normalErosion.addBands(normalRDLS)
                                  .addBands(normalWCI)
                                  .addBands(normalICHB)
                                  .addBands(normalmsdi)
                                  .addBands(normalndvi)
                                  .addBands(normalpop)
                                  .addBands(normalGDP)
var normal = combinedImage.select([0,1,2,3,4,5,6,7], normalBandNames);
// print('normalmsdi',normalmsdi)

var expressionEeq = 'normalErosion*0.03 + normalRDLS*0.06 + normalWCI*0.03 + normalICHB*0.06 + normalndvi*0.2 + normalpop*0.22 + normalGDP*0.4';

var eeq = normal.expression(expressionEeq, {
  'normalErosion': normal.select('normalErosion'),
  'normalRDLS': normal.select('normalRDLS'),
  'normalWCI': normal.select('normalWCI'),
  'normalICHB': normal.select('normalICHB'),
  'normalmsdi': normal.select('normalmsdi'),
  'normalndvi': normal.select('normalndvi'),
  'normalpop': normal.select('normalpop'),
  'normalGDP': normal.select('normalGDP'),
});


  return eeq;
}

// function MonteCarlo1(idx){
//   var rand =  Math.random()
//   rand = ee.Number(rand).float()
//   print(rand)
//   return ee.Image.constant(1).multiply(rand).clip(roi)
// }

var nq = 10

var randnum1 = ee.List([])
var randnum2 = ee.List([])
var randnum3 = ee.List([])
var randnum4 = ee.List([])
var randnum5 = ee.List([])
var randnum6 = ee.List([])
var randnum7 = ee.List([])
var randnum8 = ee.List([])
var randnum9 = ee.List([])

for (var i=0; i<nq; i++) {  
    randnum1 = randnum1.add(0.8 + Math.random() * 0.2)
    randnum2 = randnum2.add(0.8 + Math.random() * 0.2)
    randnum3 = randnum3.add(0.8 + Math.random() * 0.2)
    randnum4 = randnum4.add(0.8 + Math.random() * 0.2)
    randnum5 = randnum5.add(0.8 + Math.random() * 0.2)
    randnum6 = randnum6.add(0.8 + Math.random() * 0.2)
    randnum7 = randnum7.add(0.8 + Math.random() * 0.2)
    randnum8 = randnum8.add(0.8 + Math.random() * 0.2)
    randnum9 = randnum9.add(0.8 + Math.random() * 0.2)
}  
//print('randnum',randnum1,randnum2,randnum3,randnum4,randnum5);   

var num_list = ee.List.sequence(0,nq-1);
var imageList = num_list.map(function(idx) {return MonteCarlo(idx)});
print(imageList);
 
//myList.map(function(idx) {return ee.Image.cat([concatImage,MonteCarlo(idx, col_a1, col_lst)]);});
// var imageList = myList.map(function(idx) {return MonteCarlo1(idx);});
// print("imageList",imageList)

print(ee.ImageCollection(imageList))



//Map.addLayer(ee.ImageCollection(imageList).select(0),visRate,'test0')
//Map.addLayer(ee.ImageCollection(imageList).select(1),visRate,'test1')

var all = ee.ImageCollection(imageList)

/*



function stdCalc(image,mean) {
  var difference = image.subtract(mean)
  var stdDev = difference.multiply(difference);
  return stdDev;
}

/****************************************************************************************************************************
 * 
 *                                             Monte Carlo
 * 
 * **************************************************************************************************************************/



// //samples number
// var nq = 10;

// //var concatImage = ee.Image([]);

// var myList = ee.List.sequence(1, nq)

// //myList.map(function(idx) {return ee.Image.cat([concatImage,MonteCarlo(idx, col_a1, col_lst)]);});
// var imageList = myList.map(function(idx) {return MonteCarlo(idx);});
// print("imageList",imageList)

// Map.addLayer(ee.ImageCollection.fromImages(imageList),{},'test')


// var concatImage = ee.ImageCollection.fromImages(imageList).toBands();
// Map.addLayer(ee.ImageCollection.fromImages(imageList),{},'test')
// print("concateImage",concatImage)
 var mEEQ = all.reduce(ee.Reducer.mean());
 var stdEEQ = all.reduce(ee.Reducer.stdDev()).divide(mEEQ);


// print(stdEEQ)

// // Get the maximum value
 //var max = stdEEQ.reduceRegion({
 //  reducer: ee.Reducer.max(),
 //  geometry: roi,
 //  scale: scale_my // Specify the scale appropriate for your image
 //}); // Replace 'band_name' with the actual band name of interest

 // Print the maximum value
 //print('Maximum Value:', max);

// //print(mRSEI)
// Map.addLayer(concatImage.select(["0_constant"]),visRate,'0')
// Map.addLayer(concatImage.select(["1_constant"]),visRate,'1')
 Map.addLayer(mEEQ,visRate,'mEEQ')
 Map.addLayer(stdEEQ,visRate,'stdEEQ')

/**/
