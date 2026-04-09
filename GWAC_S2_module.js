/***************************************
 * GWAC core module for Sentinel-2
 * 
 * GWAC is a Google Earth Engine module for
 * atmospheric correction of water targets.
 *
 * For questions or feedback,
 * please contact: dianazhao0105@gmail.com
 * 
 * Versions:
 * - 2026-04-08: First public release.
 ***************************************/


/***************************************
 * Spectral constants for Sentinel-2
 ***************************************/

var bandKeys = ee.List([
  'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7',
  'B8', 'B8A', 'B9', 'B10', 'B11', 'B12'
]);

var DEFAULT_RUN_BANDS = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A'];
var S2_MAIN_BANDS = ee.List(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8']);
var S2_SWIR_BANDS = ee.List(['B8A', 'B11', 'B12']);
var S2_OUTPUT_BANDS = ee.List(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A']);
var S2_DROP_BANDS = ee.List(['B9', 'B10']);

// ----------- Band Wavelengths ----------
var A_vals = ee.List([443, 492, 560, 665, 704, 740, 783, 833, 865, 945, 1373, 1614, 2202]);
var B_vals = ee.List([442, 492, 559, 665, 704, 739, 780, 833, 864, 943, 1377, 1610, 2186]);

var bands_all = ee.Dictionary({
  'A': ee.Dictionary.fromLists(bandKeys, A_vals),
  'B': ee.Dictionary.fromLists(bandKeys, B_vals),
});
exports.bands_all = bands_all;


// ----------- Ozone ----------
A_vals = ee.List([
  0.002900, 0.026190, 0.105400, 0.050670, 0.020430,
  0.011000, 0.007081, 0.003590, 0.002180, 0.000748,
  0.000000, 0.000000, 0.000000
]);
B_vals = ee.List([
  0.002837, 0.025920, 0.104200, 0.050330, 0.020550,
  0.010960, 0.007322, 0.003559, 0.002277, 0.000746,
  0.000000, 0.000000, 0.000000
]);

var ozone_all = ee.Dictionary({
  'A': ee.Dictionary.fromLists(bandKeys, A_vals),
  'B': ee.Dictionary.fromLists(bandKeys, B_vals),
});

// ----------- Rayleigh optical depth ----------
A_vals = ee.List([
  0.235600, 0.155800, 0.090550, 0.044980, 0.035530,
  0.028970, 0.023160, 0.018530, 0.015490, 0.010830,
  0.002411, 0.001269, 0.000368
]);
B_vals = ee.List([
  0.236600, 0.156300, 0.091150, 0.044890, 0.035590,
  0.029200, 0.023530, 0.018510, 0.015540, 0.010920,
  0.002387, 0.001280, 0.000380
]);

var tau_all = ee.Dictionary({
  'A': ee.Dictionary.fromLists(bandKeys, A_vals),
  'B': ee.Dictionary.fromLists(bandKeys, B_vals),
});


/***************************************
 * S2 enrichment and angle attachment
 ***************************************/

/**
 * Enrich Sentinel-2 images with Cloud Score+ and Dynamic World layers.
 *
 * Added bands:
 * - cs_cdf
 * - label
 * - water
 *
 * Added properties:
 * - PRODUCT_URI
 * - TILE_TIME
 *
 * @param {ee.ImageCollection} S2Col
 * @param {ee.Geometry} region
 * @param {string} startDate
 * @param {string} endDate
 * @return {ee.ImageCollection}
 */
function enrichS2(S2Col, region, startDate, endDate) {
  var csPlusCol = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')
    .filterBounds(region)
    .filterDate(startDate, endDate);

  var dwCol = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
    .filterBounds(region)
    .filterDate(startDate, endDate);

  return S2Col
    .filterBounds(region)
    .filterDate(startDate, endDate)
    .map(function(img) {
      var pid = ee.String(img.get('PRODUCT_ID'));
      var parts = pid.split('_');
      return img.set({
        'PRODUCT_URI': pid.cat('.SAFE'),
        'TILE_TIME': parts.getString(5).cat('_').cat(parts.getString(2))
      });
    })
    .linkCollection(csPlusCol, ['cs_cdf'])
    .linkCollection(dwCol, ['label', 'water']);

}
exports.enrichS2 = enrichS2;


/**
 * Compute PHI from solar and view azimuth angles.
 *
 * Sentinel-2 angles are expected in degrees.
 *
 * @param {ee.Image} saa
 * @param {ee.Image} vaa
 * @return {ee.Image}
 */
function calculateRAA(saa, vaa) {
  var raa = saa.subtract(vaa).abs();
  raa = raa.where(raa.gt(180), ee.Image(360).subtract(raa));
  raa = ee.Image(180).subtract(raa);
  return raa;
}
exports.calculateRAA = calculateRAA;


/**
 * Convert HLS angular bands to GWAC-ready geometry bands.
 *
 * Added bands:
 * - PHI
 * - cosSZA
 * - cosVZA
 *
 * @param {ee.Image} hlsImg
 * @return {ee.Image}
 */
function processAngles(hlsImg) {
  var saa = hlsImg.select('SAA');
  var vaa = hlsImg.select('VAA');
  var phi = calculateRAA(saa, vaa);

  hlsImg = hlsImg
    .addBands(phi.rename('PHI'))
    .select(['SZA', 'VZA', 'PHI']);

  var cosSZA = hlsImg.select('SZA')
    .multiply(Math.PI / 180)
    .cos()
    .rename('cosSZA');

  var cosVZA = hlsImg.select('VZA')
    .multiply(Math.PI / 180)
    .cos()
    .rename('cosVZA');

  hlsImg = hlsImg.addBands(cosSZA).addBands(cosVZA);

  return hlsImg;
}
exports.processAngles = processAngles;


/**
 * Attach HLS angular information to Sentinel-2 images.
 *
 * Added bands:
 * - SZA
 * - VZA
 * - PHI
 * - cosSZA
 * - cosVZA
 *
 * Copied properties:
 * - Sensor
 * - *_scale
 *
 * @param {ee.ImageCollection} s2
 * @param {ee.Geometry} region
 * @param {string} startDate
 * @param {string} endDate
 * @return {ee.ImageCollection}
 */
function S2ToHLS(s2, region, startDate, endDate) {
  var uris = s2.aggregate_array('TILE_TIME');

  var hls = ee.ImageCollection('NASA/HLS/HLSS30/v002')
    .filterBounds(region)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.inList('system:index', uris))
    .select(['SZA', 'SAA', 'VZA', 'VAA'])
    .map(function(img) {
      var url = ee.String(img.get('PRODUCT_URI'));
      img = img.set({'Sensor': url.slice(2, 3)});
      img = processAngles(img);
      return img;
    });

  var s2HlsPairs = ee.Join.inner().apply(
    s2,
    hls,
    ee.Filter.equals({
      leftField: 'TILE_TIME',
      rightField: 'system:index'
    })
  );

  return ee.ImageCollection(s2HlsPairs.map(function(pair) {
    pair = ee.Feature(pair);

    var s2Img = ee.Image(pair.get('primary'));
    var hlsImg = ee.Image(pair.get('secondary'));

    var bBands = s2Img.bandNames()
      .filter(ee.Filter.stringStartsWith('item', 'B'));

    var extras = ee.List(['cs_cdf', 'label', 'water']);
    var existingExtras = extras.filter(
      ee.Filter.inList('item', s2Img.bandNames())
    );

    var keep = bBands.cat(existingExtras).distinct();

    var out = s2Img.select(keep).addBands(hlsImg);

    var scaleKeys = hlsImg.propertyNames()
      .filter(ee.Filter.stringStartsWith('item', 'B'));

    var fixedKeys = ee.List(['Sensor']);
    var keysToCopy = fixedKeys.cat(scaleKeys);

    return ee.Image(out.copyProperties(hlsImg, keysToCopy))
      .set({
        'image_id': s2Img.get('system:index')
      });
  }));
}
exports.S2ToHLS = S2ToHLS;


/***************************************
 * MERRA2 attachment
 ***************************************/

/**
 * Attach the closest hourly MERRA2 image to each Sentinel-2 scene.
 *
 * MERRA2 data are loaded internally from:
 * - NASA/GSFC/MERRA/slv/2
 *
 * Matching is performed within a ±1 hour time window.
 *
 * Added bands:
 * - PS
 * - TO3
 *
 * Added properties:
 * - PS_MIN
 * - PS_MAX
 *
 * @param {ee.ImageCollection} S2
 * @param {string} startDate
 * @param {string} endDate
 * @return {ee.ImageCollection}
 */
function MERRA2ToS2(S2, startDate, endDate) {
  
  var MERRA2 = ee.ImageCollection("NASA/GSFC/MERRA/slv/2")
    .filterDate(startDate, endDate)
    .select(['PS', 'TO3']);
  
  var maxDiffMs = 1 * 60 * 60 * 1000;

  var maxDiffFilter = ee.Filter.maxDifference({
    difference: maxDiffMs,
    leftField: 'system:time_start',
    rightField: 'system:time_start'
  });

  var saveBestJoin = ee.Join.saveBest({
    matchKey: 'MERRA2',
    measureKey: 'timeDiffMs'
  });

  var joined = ee.ImageCollection(
    saveBestJoin.apply(S2, MERRA2.select(['PS', 'TO3']), maxDiffFilter)
  );

  var S2_new = joined.map(function(img) {
    img = ee.Image(img);

    var geometry = img.select('B8A').geometry();
    var merra = ee.Image(img.get('MERRA2'));

    var m = merra.resample('bilinear');

    var ps = merra.select('PS');

    var psStats = ps.reduceRegion({
      reducer: ee.Reducer.minMax(),
      geometry: geometry,
      scale: 5000,
      bestEffort: true,
      maxPixels: 1e13,
      tileScale: 4
    });

    return img
      .addBands(m)
      .set({
        'PS_MIN': psStats.get('PS_min'),
        'PS_MAX': psStats.get('PS_max'),
        'MERRA2': null
      });
  });

  return S2_new;
}
exports.MERRA2ToS2 = MERRA2ToS2;


/**
 * Convenience wrapper for collection preparation.
 *
 * Workflow:
 * 1. enrichS2
 * 2. S2ToHLS
 * 3. MERRA2ToS2
 *
 * @param {ee.ImageCollection} S2Col
 * @param {ee.ImageCollection} MERRA2
 * @param {ee.Geometry} region
 * @param {string} startDate
 * @param {string} endDate
 * @return {ee.ImageCollection}
 */
function prepareS2(S2Col, region, startDate, endDate) {
  var s2 = enrichS2(S2Col, region, startDate, endDate);
  s2 = S2ToHLS(s2, region, startDate, endDate);
  s2 = MERRA2ToS2(s2, startDate, endDate);
  return s2;
}
exports.prepareS2 = prepareS2;


/***************************************
 * Water mask
 ***************************************/

/**
 * Keep clean water pixels only.
 *
 * Supported masking modes:
 * - 'ndvi' (default): cs_cdf + NDVI
 * - 'dynamicworld': cs_cdf + Dynamic World water mask
 *
 * @param {ee.Image} image
 * @param {string=} maskMethod
 * @return {ee.Image}
 */
function LandMaskS2(image, maskMethod) {
  maskMethod = maskMethod || 'ndvi';

  var keep = image.bandNames().removeAll(S2_DROP_BANDS);
  image = image.select(keep);

  var cloudSnowMask = image.select('cs_cdf').gt(0.75);

  var waterMask;
  if (maskMethod === 'dynamicworld') {
    waterMask = image.select('label').eq(0)
      .and(image.select('water').gt(0.5));
  } else {
    var ndvi = image.normalizedDifference(['B8A', 'B4']);
    waterMask = ndvi.lt(0);
  }

  var qaMask = cloudSnowMask.and(waterMask);
  return image.updateMask(qaMask);
}
exports.LandMaskS2 = LandMaskS2;


/***************************************
 * Rayleigh parameters
 ***************************************/

/**
 * Add pressure-scaled Rayleigh optical thickness and
 * per-band extraterrestrial irradiance properties.
 *
 * Required properties:
 * - Sensor
 * - SOLAR_IRRADIANCE_*
 * - REFLECTANCE_CONVERSION_CORRECTION
 *
 * @param {ee.Image} img
 * @param {ee.List} bands
 * @return {ee.Image}
 */
function getRealTau(img, bands) {
  var sp = img.select('PS');
  var sensor = img.get('Sensor');
  var tauDict = ee.Dictionary(tau_all.get(sensor));
  var dist = ee.Number(img.get('REFLECTANCE_CONVERSION_CORRECTION'));

  var out = ee.Image(ee.List(bands).iterate(function(band, acc) {
    band = ee.String(band);
    acc = ee.Image(acc);

    var tauValue = ee.Number(tauDict.get(band));
    var tauReal = ee.Image(tauValue)
      .divide(1013.25)
      .multiply(sp)
      .divide(100)
      .rename(ee.String('tau_Ray_').cat(band));

    var F0 = ee.Number(img.get(ee.String('SOLAR_IRRADIANCE_').cat(band)))
      .multiply(0.1)
      .divide(dist.pow(2));

    return acc
      .addBands(tauReal)
      .set(ee.String('F0_').cat(band), F0);
  }, img));

  return out;
}
exports.getRealTau = getRealTau;


/***************************************
 * GRAYCO model stack
 ***************************************/

// Pre-trained GRAYCO models loaded from Earth Engine assets.
var tauStart = 0;
var tauEnd = 250;
var step = 10;

var assetPrefix = 'projects/ee-dianazhao0105/assets/AC/GRAYCO_IS_';

var taus = ee.List.sequence(tauStart, tauEnd, step);

var modelList = taus.map(function(t) {
  t = ee.Number(t);
  var path = ee.String(assetPrefix).cat(t.format('%d'));
  return ee.Classifier.load(path);
});


/***************************************
 * Radiometric conversion
 ***************************************/

/**
 * Convert reflectance-like Sentinel-2 input to radiance-like
 * input for the atmospheric-correction workflow.
 *
 * Required properties:
 * - REFLECTANCE_CONVERSION_CORRECTION
 * - SOLAR_IRRADIANCE_*
 * - *_scale
 *
 * @param {ee.Image} img
 * @param {ee.List} bands
 * @return {ee.Image}
 */
function rhoToLt(img, bands) {
  var cosSZA = img.select('cosSZA');

  var imgs = ee.Image(ee.List(bands).iterate(function(band, acc) {
    band = ee.String(band);
    acc = ee.Image(acc);

    var dist = ee.Number(img.get('REFLECTANCE_CONVERSION_CORRECTION'));
    var scale = ee.Number(img.get(ee.String('SOLAR_IRRADIANCE_').cat(band)))
      .multiply(0.1)
      .divide(dist.pow(2));

    scale = scale
      .multiply(ee.Number(img.get(band.cat('_scale'))))
      .divide(Math.PI);

    var one = img.select([band])
      .multiply(cosSZA)
      .multiply(scale)
      .rename([band]);

    return acc.addBands(one);
  }, ee.Image([])));

  return img.addBands(imgs, null, true);
}
exports.rhoToLt = rhoToLt;


/***************************************
 * Atmospheric transmittance
 ***************************************/

function transOzone(tauOzone, cosAngle) {
  var tau = cosAngle.multiply(0).add(tauOzone);
  return tau.multiply(-1).divide(cosAngle).exp();
}

function transRay(tauRay, cosAngle) {
  var tau = cosAngle.multiply(0).add(tauRay);
  return tau.multiply(-0.5).divide(cosAngle).exp();
}


/***************************************
 * Gas correction
 ***************************************/

/**
 * Apply ozone correction to the selected Sentinel-2 bands.
 *
 * @param {ee.Image} img
 * @param {ee.List} bands
 * @return {ee.Image}
 */
function gasCorrection(img, bands) {
  var cosSZA = img.select('cosSZA');
  var cosVZA = img.select('cosVZA');

  var sensor = img.get('Sensor');
  var ozoneDict = ee.Dictionary(ozone_all.get(sensor));
  var ozoneScale = img.select('TO3').multiply(0.001);

  var imgs = ee.Image(ee.List(bands).iterate(function(band, acc) {
    band = ee.String(band);
    acc = ee.Image(acc);

    var ozoneBand = ozoneScale.multiply(ee.Number(ozoneDict.get(band)));

    var tOzone = transOzone(ozoneBand, cosSZA);
    tOzone = tOzone.multiply(transOzone(ozoneBand, cosVZA));

    var one = img.select([band]).divide(tOzone).rename([band]);
    return acc.addBands(one);
  }, ee.Image([])));

  return img.addBands(imgs, null, true);
}
exports.gasCorrection = gasCorrection;


/***************************************
 * Rayleigh scattering with GRAYCO
 ***************************************/

/**
 * Compute the angular factor used by GRAYCO.
 *
 * Input angle units: degree × 10000.
 *
 * @param {ee.Image} angles
 * @return {ee.Image}
 */
function rayleighFactor(angles) {
  var t0x = angles.select('S');
  var tvx = angles.select('V');
  var dpx = angles.select('P');

  var deg100_to_rad = ee.Image.constant(Math.PI / 180.0 / 10000.0);
  var EPS = ee.Image.constant(1e-6);
  var ONE = ee.Image.constant(1.0);
  var C075 = ee.Image.constant(0.75);
  var FOUR = ee.Image.constant(4.0);
  var PI = ee.Image.constant(Math.PI);

  var t0 = t0x.multiply(deg100_to_rad);
  var tv = tvx.multiply(deg100_to_rad);
  var dp = dpx.multiply(deg100_to_rad);

  var c0 = t0.cos();
  var s0 = t0.sin();
  var cv = tv.cos();
  var sv = tv.sin();

  var cdp = dp.cos().multiply(-1);

  var c0cv = c0.multiply(cv);
  var s0sv = s0.multiply(sv);

  var ct_plus = c0cv.subtract(s0sv.multiply(cdp));
  var ct_minus = c0cv.multiply(-1).subtract(s0sv.multiply(cdp));

  var Pplus = C075.multiply(ONE.add(ct_plus.multiply(ct_plus)));
  var Pminus = C075.multiply(ONE.add(ct_minus.multiply(ct_minus)));

  var nw = ee.Number(1.34);

  function fresnel_r(theta) {
    var sinT = theta.sin();
    var tt = sinT.multiply(nw).clamp(-1, 1).asin();

    var tPlus = theta.add(tt);
    var tMinus = theta.subtract(tt);

    var sinP = tPlus.sin();
    var sinM = tMinus.sin();
    var cosP = tPlus.cos();
    var cosM = tMinus.cos();

    var rs = sinM.divide(sinP.add(EPS));
    rs = rs.multiply(rs);

    var rp = sinM.multiply(cosP)
      .divide(cosM.multiply(sinP).add(EPS));
    rp = rp.multiply(rs);

    return rs.add(rp).multiply(0.5);
  }

  var rsum = fresnel_r(tv).add(fresnel_r(t0));
  var pr = Pminus.add(rsum.multiply(Pplus));

  var denom = FOUR.multiply(cv.abs().add(EPS)).multiply(PI);
  return pr.divide(denom);
}


/**
 * Remove Rayleigh path radiance for visible/NIR bands
 * using pre-trained GRAYCO models.
 *
 * @param {ee.Image} img
 * @param {ee.List} bands
 * @return {ee.Image}
 */
function getLrc(img, bands) {
  var one = ee.Image.constant(1);
  var waterMask = img.select('B8A').mask();

  var bandInput = img.select(['SZA', 'VZA', 'PHI'])
    .rename(['S', 'V', 'P'])
    .updateMask(waterMask)
    .multiply(10000);

  var factor = rayleighFactor(bandInput);
  var I_S = factor.rename('I_S');

  var tauScale = ee.Number(100);
  var scale1000 = tauScale.divide(1000);
  var kStep = ee.Number(1).divide(scale1000).int();

  var psMin = ee.Number(ee.Algorithms.If(img.get('PS_MIN'), img.get('PS_MIN'), 101325));
  var psMax = ee.Number(ee.Algorithms.If(img.get('PS_MAX'), img.get('PS_MAX'), 101325));

  psMin = psMin.divide(100).divide(1013.25);
  psMax = psMax.divide(100).divide(1013.25);

  var sensor = img.get('Sensor');
  var tauDict = ee.Dictionary(tau_all.get(sensor));

  var bandNames = S2_MAIN_BANDS.filter(ee.Filter.inList('item', ee.List(bands)));

  var bandResults = bandNames.map(function(band) {
    band = ee.String(band);
    var bandName = ee.String('tau_Ray_').cat(band);

    var bandTau = img.select(bandName).multiply(1000);

    var tauIdx = bandTau.multiply(scale1000).floor()
      .multiply(1000).divide(tauScale).int();
    tauIdx = tauIdx.clamp(0, 250);
    bandTau = bandTau.clamp(0, 250);

    var tauConst = ee.Number(tauDict.get(band));

    var tau1000Min = tauConst.multiply(psMin).multiply(1000).clamp(0, 250);
    var tau1000Max = tauConst.multiply(psMax).multiply(1000).clamp(0, 250);

    var tau1000ToIdx = function(t) {
      t = ee.Number(t);
      return t.multiply(scale1000).floor()
        .multiply(1000).divide(tauScale).int()
        .clamp(0, 250);
    };

    var kMin = tau1000ToIdx(tau1000Min);
    var kMax = tau1000ToIdx(tau1000Max);

    var kkMin = kMin.min(kMax);
    var kkMax = kMin.max(kMax);

    var uniqueTau = ee.List.sequence(kkMin, kkMax, 10);
    var acc0 = img.select('SZA').multiply(0).rename('acc');

    var accumulated = ee.Image(uniqueTau.iterate(function(k, accImg) {
      k = ee.Number(k);
      var kNext = k.add(kStep).min(250);

      var mask = tauIdx.eq(k).selfMask();
      accImg = ee.Image(accImg);

      var interInput = bandInput.updateMask(mask);
      interInput = interInput.addBands(I_S.updateMask(mask));

      var idx1 = ee.Number(k).divide(10).int();
      var idx2 = ee.Number(kNext).divide(10).int();

      var model_1 = modelList.get(idx1);
      var model_2 = modelList.get(idx2);

      var res1 = interInput.classify(model_1).toFloat();
      var res2 = interInput.classify(model_2).toFloat();

      var coef1 = bandTau.updateMask(mask).subtract(ee.Image(k));
      coef1 = coef1.multiply(scale1000);

      var blended = res1.multiply(one.subtract(coef1))
        .add(res2.multiply(coef1));

      return accImg.where(mask, blended);
    }, acc0));

    var F0Band = ee.Number(img.get(ee.String('F0_').cat(band)));
    var resultPerBand = accumulated.multiply(F0Band);
    resultPerBand = img.select([band]).subtract(resultPerBand);

    return ee.Image(resultPerBand).rename([band]);
  });

  bandResults = ee.Image(bandResults.iterate(function(bandImg, acc) {
    return ee.Image(acc).addBands(ee.Image(bandImg));
  }, ee.Image([])));

  return img.addBands(bandResults, null, true);
}
exports.getLrc = getLrc;


/**
 * Remove Rayleigh path radiance for B8A, B11, and B12
 * using pre-trained GRAYCO models.
 *
 * @param {ee.Image} img
 * @param {ee.List} bands
 * @return {ee.Image}
 */
function SWIRLrc(img, bands) {
  var one = ee.Image.constant(1);
  var waterMask = img.select('B8A').mask();

  var bandInput = img.select(['SZA', 'VZA', 'PHI'])
    .rename(['S', 'V', 'P'])
    .updateMask(waterMask)
    .multiply(10000);

  var factor = rayleighFactor(bandInput);
  var I_S = factor.rename('I_S');

  var tauScale = ee.Number(100);
  var scale1000 = tauScale.divide(1000);
  var kStep = ee.Number(1).divide(scale1000).int();

  var psMin = ee.Number(ee.Algorithms.If(img.get('PS_MIN'), img.get('PS_MIN'), 101325));
  var psMax = ee.Number(ee.Algorithms.If(img.get('PS_MAX'), img.get('PS_MAX'), 101325));

  psMin = psMin.divide(100).divide(1013.25);
  psMax = psMax.divide(100).divide(1013.25);

  var sensor = img.get('Sensor');
  var tauDict = ee.Dictionary(tau_all.get(sensor));

  var bandNames = S2_SWIR_BANDS.filter(ee.Filter.inList('item', ee.List(bands)));

  var bandResults = bandNames.map(function(band) {
    band = ee.String(band);
    var bandName = ee.String('tau_Ray_').cat(band);

    var bandTau = img.select(bandName).multiply(1000);

    var tauIdx = bandTau.multiply(scale1000).floor()
      .multiply(1000).divide(tauScale).int();
    tauIdx = tauIdx.clamp(0, 250);
    bandTau = bandTau.clamp(0, 250);

    var tauConst = ee.Number(tauDict.get(band));

    var tau1000Min = tauConst.multiply(psMin).multiply(1000).clamp(0, 250);
    var tau1000Max = tauConst.multiply(psMax).multiply(1000).clamp(0, 250);

    var tau1000ToIdx = function(t) {
      t = ee.Number(t);
      return t.multiply(scale1000).floor()
        .multiply(1000).divide(tauScale).int()
        .clamp(0, 250);
    };

    var kMin = tau1000ToIdx(tau1000Min);
    var kMax = tau1000ToIdx(tau1000Max);

    var kkMin = kMin.min(kMax);
    var kkMax = kMin.max(kMax);

    var uniqueTau = ee.List.sequence(kkMin, kkMax, 10);
    var acc0 = img.select('SZA').multiply(0).rename('acc');

    var accumulated = ee.Image(uniqueTau.iterate(function(k, accImg) {
      k = ee.Number(k);
      var kNext = k.add(kStep).min(250);

      var mask = tauIdx.eq(k).selfMask();
      accImg = ee.Image(accImg);

      var interInput = bandInput.updateMask(mask);
      interInput = interInput.addBands(I_S.updateMask(mask));

      var idx1 = ee.Number(k).divide(10).int();
      var idx2 = ee.Number(kNext).divide(10).int();

      var model_1 = modelList.get(idx1);
      var model_2 = modelList.get(idx2);

      var res1 = interInput.classify(model_1).toFloat();
      var res2 = interInput.classify(model_2).toFloat();

      var coef1 = bandTau.updateMask(mask).subtract(ee.Image(k));
      coef1 = coef1.multiply(scale1000);

      var blended = res1.multiply(one.subtract(coef1))
        .add(res2.multiply(coef1));

      return accImg.where(mask, blended);
    }, acc0));

    var F0Band = ee.Number(img.get(ee.String('F0_').cat(band)));
    var resultPerBand = accumulated.multiply(F0Band);
    resultPerBand = img.select([band]).subtract(resultPerBand);

    return ee.Image(resultPerBand).rename([band]);
  });

  bandResults = ee.Image(bandResults.iterate(function(bandImg, acc) {
    return ee.Image(acc).addBands(ee.Image(bandImg));
  }, ee.Image([])));

  return img.addBands(bandResults, null, true);
}
exports.SWIRLrc = SWIRLrc;


/***************************************
 * SWIR extrapolation
 ***************************************/

/**
 * Estimate the Angstrom-like slope term used in the SWIR extrapolation step.
 *
 * @param {ee.Image} img
 * @return {ee.Image}
 */
function getAngstrom(img) {
  var sensor = img.get('Sensor');
  var wlDict = ee.Dictionary(bands_all.get(sensor));

  var wl8A = ee.Number(wlDict.get('B8A'));
  var wl11 = ee.Number(wlDict.get('B11'));
  var wl12 = ee.Number(wlDict.get('B12'));

  var f0_8A = ee.Image.constant(ee.Number(img.get('SOLAR_IRRADIANCE_B8A')));
  var f0_11 = ee.Image.constant(ee.Number(img.get('SOLAR_IRRADIANCE_B11')));
  var f0_12 = ee.Image.constant(ee.Number(img.get('SOLAR_IRRADIANCE_B12')));

  var eps = 1e-10;

  var n8A = img.select('B8A').divide(f0_8A).max(eps);
  var n11 = img.select('B11').divide(f0_11).max(eps);
  var n12 = img.select('B12').divide(f0_12).max(eps);

  var log11_12 = n11.divide(n12).log();
  var log8A_12 = n8A.divide(n12).log();

  var slope11 = log11_12.divide(wl12.subtract(wl11));
  var slope8A = log8A_12.divide(wl12.subtract(wl8A));
  var am_strong = slope11.min(slope8A).rename('am_strong');

  var all = img.bandNames();
  all = all.removeAll([
    'B9', 'B10', 'B11',
    'tau_Ray_B11', 'tau_Ray_B12',
    'Lr_B11', 'Lr_B12'
  ]);

  return img.select(all).addBands(am_strong);
}
exports.getAngstrom = getAngstrom;


/**
 * Extrapolate aerosol path radiance from B12 and remove it
 * from the output bands.
 *
 * @param {ee.Image} img
 * @param {ee.List} runBands
 * @return {ee.Image}
 */
function getLw(img, runBands) {
  var sensor = img.get('Sensor');
  var wlDict = ee.Dictionary(bands_all.get(sensor));

  var wlB12 = ee.Number(wlDict.get('B12'));
  var f0B12 = ee.Number(img.get('SOLAR_IRRADIANCE_B12'));
  var LaB12 = img.select(['B12']);

  var imgs = ee.Image(ee.List(runBands).iterate(function(band, acc) {
    band = ee.String(band);
    acc = ee.Image(acc);

    var wlBand = ee.Number(wlDict.get(band));
    var f0Band = ee.Number(img.get(ee.String('SOLAR_IRRADIANCE_').cat(band)));

    var La = LaB12.multiply(f0Band).divide(f0B12)
      .multiply(
        img.select(['am_strong'])
          .multiply(wlB12.subtract(wlBand))
          .exp()
      );

    var one = img.select([band]).subtract(La).rename([band]);
    return acc.addBands(one);
  }, ee.Image([])));

  return img.addBands(imgs, null, true);
}
exports.getLw = getLw;


/***************************************
 * Final Rrs
 ***************************************/

/**
 * Compute remote-sensing reflectance for the output bands.
 *
 * Returned band names follow the original Sentinel-2 band names.
 *
 * @param {ee.Image} img
 * @param {ee.List} runBands
 * @return {ee.Image}
 */
function getRrs(img, runBands) {
  var cosSZA = img.select('cosSZA');
  var cosVZA = img.select('cosVZA');

  var imgs = ee.Image(ee.List(runBands).iterate(function(band, acc) {
    band = ee.String(band);
    acc = ee.Image(acc);

    var F0Band = ee.Number(img.get(ee.String('F0_').cat(band)));
    var Ray = img.select([ee.String('tau_Ray_').cat(band)]);

    var tsRay = transRay(Ray, cosSZA);
    var tvRay = transRay(Ray, cosVZA);

    var one = img.select([band])
      .divide(tsRay)
      .divide(tvRay)
      .divide(cosSZA)
      .divide(F0Band)
      .rename([band]);

    return acc.addBands(one);
  }, ee.Image([])));

  return ee.Image(imgs.copyProperties(img));
}
exports.getRrs = getRrs;


/***************************************
 * Band selection
 ***************************************/

/**
 * Build the Sentinel-2 processing band groups.
 *
 * Output bands:
 * - requested run bands, excluding B9, B10, B11, and B12
 *
 * Processing bands:
 * - output bands + B8A + B11 + B12
 *
 * @param {Array|ee.List=} runBands
 * @return {ee.Dictionary}
 */
function buildBandsToProcess(runBands) {
  runBands = runBands || DEFAULT_RUN_BANDS;
  runBands = ee.List(runBands).removeAll(['B9', 'B10', 'B11', 'B12']);

  var must = ee.List(['B8A', 'B11', 'B12']);
  var bandsToProcess = runBands.cat(must).distinct();

  return ee.Dictionary({
    run: runBands,
    all: bandsToProcess
  });
}
exports.buildBandsToProcess = buildBandsToProcess;


/***************************************
 * Main workflow
 ***************************************/

/**
 * Run the GWAC workflow on a prepared Sentinel-2 image.
 *
 * Required preprocessing:
 * 1. enrichS2()
 * 2. S2ToHLS()
 * 3. MERRA2ToS2()
 * or JUST prepareS2()
 *
 * Default output bands:
 * - B1-B8A, excluding B9, B10, B11, and B12
 *
 * Processing bands:
 * - output bands + B8A + B11 + B12
 *
 * Supported masking modes:
 * - 'ndvi' (default)
 * - 'dynamicworld'
 *
 * @param {ee.Image} img
 * @param {Array|ee.List|string=} runBands
 * @param {string=} maskMethod
 * @return {ee.Image}
 */
function GWAC(img, runBands, maskMethod) {
  if (typeof runBands === 'string' && maskMethod === undefined) {
    maskMethod = runBands;
    runBands = null;
  }

  maskMethod = maskMethod || 'ndvi';

  var bd = buildBandsToProcess(runBands);
  var outBands = ee.List(bd.get('run'));
  var procBands = ee.List(bd.get('all'));

  img = LandMaskS2(img, maskMethod);
  img = getRealTau(img, procBands);
  img = rhoToLt(img, procBands);
  img = gasCorrection(img, procBands);

  img = getLrc(img, procBands);
  img = SWIRLrc(img, procBands);
  img = getAngstrom(img);
  img = getLw(img, outBands);

  return getRrs(img, outBands);
}
exports.GWAC = GWAC;