/***************************************
 * GWAC core module for Landsat 8
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
 * Spectral constants for Landsat 8
 ***************************************/

var bandKeys = ee.List(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7']);
var bandWaves = ee.List([443, 482, 561, 665, 865, 1609, 2201]);

var bands_all = ee.Dictionary.fromLists(bandKeys, bandWaves);
exports.bands_all = bands_all;

// Rayleigh Optical Depth
var tau = ee.List([
  2.352E-01, 1.685E-01, 9.021E-02, 4.794E-02,
  1.551E-02, 1.284E-03, 3.697E-04
]);
var tau_all = ee.Dictionary.fromLists(bandKeys, tau);

// Ozone
var k_oz = ee.List([
  2.929E-03, 1.957E-02, 1.038E-01, 6.201E-02,
  2.223E-03, 0.000E+00, 0.000E+00
]);
var ozone_all = ee.Dictionary.fromLists(bandKeys, k_oz);

// F0
var F0 = ee.List([
  1894.746, 2004.869, 1820.121, 1550.348,
  951.206, 247.557, 85.460
]);
var F0_all = ee.Dictionary.fromLists(bandKeys, F0);


/***************************************
 * MERRA2 attachment
 ***************************************/

/**
 * Attach the closest hourly MERRA2 image to each Landsat scene.
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
 * @param {ee.ImageCollection} Landsat
 * @param {string} startDate
 * @param {string} endDate
 * @return {ee.ImageCollection}
 */
function MERRA2ToL89(L89, startDate, endDate) {
  
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
    saveBestJoin.apply(L89, MERRA2.select(['PS', 'TO3']), maxDiffFilter)
  );

  var L89_new = joined.map(function(img) {
    img = ee.Image(img);

    var geometry = img.select('B2').geometry();
    var merra = ee.Image(img.get('MERRA2'));

    var m = merra
      .clip(geometry.buffer(150000))
      .resample('bilinear')
      .clip(geometry);

    var ps = merra.select('PS');

    var psStats = ps.reduceRegion({
      reducer: ee.Reducer.minMax(),
      geometry: geometry,
      scale: 5000,
      bestEffort: true,
      maxPixels: 1e13,
      tileScale: 4
    });

    var psMin = psStats.get('PS_min');
    var psMax = psStats.get('PS_max');

    return img
      .addBands(m)
      .set({'PS_MIN': psMin, 'PS_MAX': psMax})
      .set('MERRA2', null);
  });

  return L89_new;
}
exports.MERRA2ToL89 = MERRA2ToL89;


/***************************************
 * Water mask
 ***************************************/

/**
 * Keep clean water pixels only.
 *
 * Rules:
 * - QA_PIXEL bit 0-5 must all be 0
 * - QA_PIXEL bit 7 must be 1 (water)
 * - QA_RADSAT must be 0
 *
 * @param {ee.Image} image
 * @return {ee.Image}
 */
function LandMaskL89(image) {
  var drop = ee.List(['B8', 'B9', 'B10', 'B11']);
  var keep = image.bandNames().removeAll(drop);
  image = image.select(keep);

  var qa = image.select('QA_PIXEL');
  var satMask = image.select('QA_RADSAT').eq(0);

  var cloudSnowMask = qa.bitwiseAnd(63).eq(0);
  var waterMask = qa.bitwiseAnd(1 << 7).gt(0);

  var qaMask = cloudSnowMask.and(waterMask).and(satMask);
  var qaKeep = image.select(['QA_PIXEL', 'QA_RADSAT']);

  image = image.updateMask(qaMask);

  return image.addBands(qaKeep, null, true);
}
exports.LandMaskL89 = LandMaskL89;


/***************************************
 * Geometry helpers
 ***************************************/

/**
 * Compute PHI from solar and view azimuth angles.
 *
 * Landsat 8 angles are stored as degree × 100.
 *
 * @param {ee.Image} saa
 * @param {ee.Image} vaa
 * @return {ee.Image}
 */
function calculateRAA(saa, vaa) {
  var raa = saa.subtract(vaa).abs();

  raa = raa.where(raa.gt(18000), ee.Image(36000).subtract(raa));
  raa = ee.Image(18000).subtract(raa);

  return raa;
}


/***************************************
 * Rayleigh parameters
 ***************************************/

/**
 * Add pressure-scaled Rayleigh optical thickness and
 * per-band extraterrestrial irradiance properties.
 *
 * @param {ee.Image} img
 * @param {ee.List} bands
 * @return {ee.Image}
 */
function getRealTau(img, bands) {
  var saa = img.select('SAA');
  var vaa = img.select('VAA');
  var raa = calculateRAA(saa, vaa);

  img = img.select(img.bandNames().removeAll(['SAA', 'VAA']))
    .addBands(raa.rename('PHI'));

  var sp = img.select('PS');
  var dist = ee.Number(img.get('EARTH_SUN_DISTANCE'));

  var out = ee.Image(bands.iterate(function(band, acc) {
    band = ee.String(band);
    acc = ee.Image(acc);

    var bName = ee.String('tau_Ray_').cat(band);
    var tau_value = ee.Number(tau_all.get(band));

    var tau_real = sp.divide(100)
      .divide(1013.25)
      .multiply(tau_value)
      .rename(bName);

    var F0_band = ee.Number(F0_all.get(band))
      .multiply(0.1)
      .divide(dist.pow(2));

    return acc
      .addBands(tau_real)
      .set(ee.String('F0_').cat(band), F0_band);

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
 * Convert DN to TOA spectral radiance for Landsat 8 C2 L1.
 *
 * @param {ee.Image} img
 * @param {ee.List} bands
 * @return {ee.Image}
 */
function DNToLt(img, bands) {
  img = ee.Image(img);

  var Lt = ee.Image(bands.iterate(function(band, acc) {
    band = ee.String(band);
    acc = ee.Image(acc);

    var nStr = band.slice(1);

    var ML = ee.Number(img.get(ee.String('RADIANCE_MULT_BAND_').cat(nStr)));
    var AL = ee.Number(img.get(ee.String('RADIANCE_ADD_BAND_').cat(nStr)));

    var lt = img.select(band).multiply(ML).add(AL)
      .multiply(0.1);

    return acc.addBands(lt);
  }, ee.Image([])));

  return img.addBands(Lt, null, true);
}
exports.DNToLt = DNToLt;


/***************************************
 * Atmospheric transmittance
 ***************************************/

function transOzone(tauOzone, cosAngle) {
  var tauImg = cosAngle.multiply(0).add(tauOzone);
  return tauImg.multiply(-1).divide(cosAngle).exp();
}

function transRay(tauRay, cosAngle) {
  var tauImg = cosAngle.multiply(0).add(tauRay);
  return tauImg.multiply(-0.5).divide(cosAngle).exp();
}


/***************************************
 * Gas correction
 ***************************************/

/**
 * Apply ozone correction to the selected bands.
 *
 * @param {ee.Image} img
 * @param {ee.List} bands
 * @return {ee.Image}
 */
function gasCorrection(img, bands) {
  var cosSZA = img.select('SZA').divide(100).multiply(Math.PI / 180).cos();
  var cosVZA = img.select('VZA').divide(100).multiply(Math.PI / 180).cos();

  var ozoneScale = img.select('TO3').multiply(0.001);

  var imgs = ee.Image(bands.iterate(function(band, acc) {
    band = ee.String(band);
    acc = ee.Image(acc);

    var ozoneBand = ozoneScale.multiply(ee.Number(ozone_all.get(band)));

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
  var factor = pr.divide(denom);

  return factor;
}


/**
 * Remove Rayleigh path radiance using pre-trained GRAYCO models.
 *
 * @param {ee.Image} img
 * @param {ee.List} bands
 * @return {ee.Image}
 */
function getLrc(img, bands) {
  var one = ee.Image.constant(1);
  var waterMask = img.select('B2').mask();

  var bandInput = img.select(['SZA', 'VZA', 'PHI'])
    .rename(['S', 'V', 'P'])
    .updateMask(waterMask)
    .multiply(100);

  var factor = rayleighFactor(bandInput);
  var I_S = factor.rename('I_S');

  var tauScale = ee.Number(100);
  var scale1000 = tauScale.divide(1000);
  var kStep = ee.Number(1).divide(scale1000).int();

  var psMin = ee.Number(ee.Algorithms.If(img.get('PS_MIN'), img.get('PS_MIN'), 101325));
  var psMax = ee.Number(ee.Algorithms.If(img.get('PS_MAX'), img.get('PS_MAX'), 101325));

  psMin = psMin.divide(100).divide(1013.25);
  psMax = psMax.divide(100).divide(1013.25);

  var bandResults = bands.map(function(band) {
    var bandName = ee.String('tau_Ray_').cat(ee.String(band));

    var bandTau = img.select(bandName).multiply(1000);

    var tauIdx = bandTau.multiply(scale1000).floor()
      .multiply(1000).divide(tauScale).int();
    tauIdx = tauIdx.clamp(0, 250);
    bandTau = bandTau.clamp(0, 250);

    var tauConst = ee.Number(tau_all.get(ee.String(band)));

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

    var acc0 = img.select('B2').multiply(0).rename('acc');

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

      var accWhenHasPix = accImg.where(mask, blended);

      return accWhenHasPix;
    }, acc0));

    var dist = ee.Number(img.get(ee.String('F0_').cat(band)));
    var resultPerBand = accumulated.multiply(dist);

    resultPerBand = img.select([band]).subtract(resultPerBand);

    return ee.Image(resultPerBand).rename([band]);
  });

  bandResults = ee.Image(bandResults.iterate(function(bandImg, acc) {
    return ee.Image(acc).addBands(ee.Image(bandImg));
  }, ee.Image([])));

  img = img.addBands(bandResults, null, true);

  return img;
}
exports.getLrc = getLrc;


/***************************************
 * SWIR extrapolation
 ***************************************/

var wl5 = ee.Number(bands_all.get('B5'));
var wl6 = ee.Number(bands_all.get('B6'));
var wl7 = ee.Number(bands_all.get('B7'));

var f0_5 = ee.Image.constant(ee.Number(F0_all.get('B5')));
var f0_6 = ee.Image.constant(ee.Number(F0_all.get('B6')));
var f0_7 = ee.Image.constant(ee.Number(F0_all.get('B7')));

/**
 * Estimate the Angstrom-like slope term used in the SWIR extrapolation step.
 *
 * @param {ee.Image} img
 * @return {ee.Image}
 */
function getAngstrom(img) {
  var n5 = img.select('B5').divide(f0_5);
  var n6 = img.select('B6').divide(f0_6);
  var n7 = img.select('B7').divide(f0_7);

  var log6_7 = n6.divide(n7).log();
  var log5_7 = n5.divide(n7).log();

  var slope6 = log6_7.divide(wl7.subtract(wl6));
  var slope5 = log5_7.divide(wl7.subtract(wl5));
  var am_strong = slope6.min(slope5).rename('am_strong');

  var all = img.bandNames();
  all = all.removeAll(['B6', 'tau_Ray_B6']);

  return img.select(all).addBands(am_strong);
}
exports.getAngstrom = getAngstrom;


/**
 * Extrapolate aerosol path radiance from B7 and remove it
 * from the requested output bands.
 *
 * @param {ee.Image} img
 * @param {ee.List} runBands
 * @return {ee.Image}
 */
function getLw(img, runBands) {
  var f0B7 = ee.Number(img.get(ee.String('F0_B7')));
  var LaB7 = img.select(['B7']);

  var imgs = ee.Image(runBands.iterate(function(band, acc) {
    band = ee.String(band);
    acc = ee.Image(acc);

    var wlBand = ee.Number(bands_all.get(band));
    var f0Band = ee.Number(img.get(ee.String('F0_').cat(band)));

    var La = LaB7.multiply(f0Band).divide(f0B7)
      .multiply(img.select(['am_strong']).multiply(
        wl7.subtract(wlBand)).exp());

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
 * Compute remote-sensing reflectance for the requested bands.
 *
 * @param {ee.Image} img
 * @param {ee.List} runBands
 * @return {ee.Image}
 */
function getRrs(img, runBands) {
  var cosSZA = img.select('SZA').divide(100).multiply(Math.PI / 180).cos();
  var cosVZA = img.select('VZA').divide(100).multiply(Math.PI / 180).cos();

  var imgs = ee.Image(runBands.iterate(function(band, acc) {
    band = ee.String(band);
    acc = ee.Image(acc);

    var F0Band = ee.Number(img.get(ee.String('F0_').cat(band)));

    var Ray = img.select([ee.String('tau_Ray_').cat(band)]);
    var tsRay = transRay(Ray, cosSZA);
    var tvRay = transRay(Ray, cosVZA);

    var one = img.select([band]).divide(tsRay).divide(tvRay)
      .divide(cosSZA).divide(F0Band).rename([band]);

    return acc.addBands(one);
  }, ee.Image([])));

  return imgs;
}
exports.getRrs = getRrs;


/***************************************
 * Band selection
 ***************************************/

/**
 * Build the processing band list.
 *
 * Output bands:
 * - requested run bands, excluding B6 and B7
 *
 * Processing bands:
 * - output bands + B5 + B6 + B7
 *
 * @param {Array|ee.List=} runBands
 * @return {ee.Dictionary}
 */
function buildBandsToProcess(runBands) {
  runBands = runBands || ['B1', 'B2', 'B3', 'B4', 'B5'];
  runBands = ee.List(runBands).removeAll(['B6', 'B7']);

  var must = ee.List(['B5', 'B6', 'B7']);
  var bandsToProcess = runBands.cat(must).distinct();

  return ee.Dictionary({
    run: runBands,
    all: bandsToProcess
  });
}


/***************************************
 * Main workflow
 ***************************************/

/**
 * Run the GWAC workflow on Landsat 8 input.
 *
 * Required preprocessing:
 * MERRA2ToL89()
 *
 * @param {ee.Image} img
 * @param {Array|ee.List=} runBands
 * @return {ee.Image}
 */
function acSWIR(img, runBands) {
  runBands = runBands || ['B1', 'B2', 'B3', 'B4', 'B5'];

  var bd = buildBandsToProcess(runBands);

  var outBands = ee.List(bd.get('run'));
  var procBands = ee.List(bd.get('all'));

  img = LandMaskL89(img);
  img = getRealTau(img, procBands);
  img = DNToLt(img, procBands);
  img = gasCorrection(img, procBands);

  img = getLrc(img, procBands);
  img = getAngstrom(img);
  img = getLw(img, outBands);
  img = getRrs(img, outBands);

  return img;
}
exports.GWAC = acSWIR;