const Point = (function(){
  function P(x, y){
    this.x = x;
    this.y = y;
  }
  return P;
})();

const LatLon = (function(){
  const p = "prototype",
    m = Math,
    mcos = m.cos,
    msin = m.sin,
    mpow = m.pow,
    msqr = m.sqrt;
  function L(lat, lon){
    this.lat = lat;
    this.lon = lon;
  }
  L[p].distance = function(latlon){
    const C = m.PI/180,
      dlat = this.lat - latlon.lat,
      dlon = this.lon - latlon.lon,
      a = mpow(msin(C*dlat / 2), 2) + mcos(C*this.lat) * mcos(C*latlon.lat) * mpow(msin(C*dlon / 2), 2);
    return 12756274 * m.atan2(msqr(a), msqr(1 - a));
  };
  return L;
})();

const SpheroidProjection = (function(){
  const p = "prototype",
    m = Math,
    pi = m.PI,
    _180 = 180.0,
    rad = 6378137,
    originShift = pi * rad,
    pi_180 = pi/_180;
  function S(){
  }
  S[p].latLonToMeters = function(latlon){
    return new Point(
      latlon.lon*rad*pi_180,
       m.log(m.tan((90+latlon.lat)*pi_180/2))*rad
    );
  };
  S[p].metersToLatLon = function(mxy){
    return new LatLon(
      (2*m.atan(m.exp(mxy.y/rad))-pi/2)/pi_180,
       mxy.x/rad/pi_180
    );
  };
  S[p].resolution = function(zoom){
    return (2 * originShift) / (256 * m.pow(2, zoom));
  };
  S[p].zoomForPixelSize = function(pixelSize ){
    for(let i=0; i<30; i++){
      if(pixelSize > this.resolution(i)){
        return m.max(i-1,0);
      }
    }
  };
  S[p].pixelsToMeters = function(px, py, zoom){
    const res = this.resolution( zoom ),
      mx = px * res - originShift,
      my = py * res - originShift;
    return new Point(mx, my);
  };
  return S;
})();


function adj(m) { // Compute the adjugate of m
    return [
        m[4] * m[8] - m[5] * m[7], m[2] * m[7] - m[1] * m[8], m[1] * m[5] - m[2] * m[4],
        m[5] * m[6] - m[3] * m[8], m[0] * m[8] - m[2] * m[6], m[2] * m[3] - m[0] * m[5],
        m[3] * m[7] - m[4] * m[6], m[1] * m[6] - m[0] * m[7], m[0] * m[4] - m[1] * m[3]
    ];
}
function multmm(a, b) { // multiply two matrices
    var c = Array(9);
    for (var i = 0; i !== 3; ++i) {
        for (var j = 0; j !== 3; ++j) {
            var cij = 0;
            for (var k = 0; k !== 3; ++k) {
                cij += a[3 * i + k] * b[3 * k + j];
            }
            c[3 * i + j] = cij;
        }
    }
    return c;
}
function multmv(m, v) { // multiply matrix and vector
    return [
        m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
        m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
        m[6] * v[0] + m[7] * v[1] + m[8] * v[2]
    ];
}
function basisToPoints(x1, y1, x2, y2, x3, y3, x4, y4) {
    var m = [
        x1, x2, x3,
        y1, y2, y3,
        1,  1,  1
    ];
    var v = multmv(adj(m), [x4, y4, 1]);
    return multmm(m, [
        v[0], 0, 0,
        0, v[1], 0,
        0, 0, v[2]
    ]);
}
function general2DProjection(
    x1s, y1s, x1d, y1d,
    x2s, y2s, x2d, y2d,
    x3s, y3s, x3d, y3d,
    x4s, y4s, x4d, y4d
) {
    var s = basisToPoints(x1s, y1s, x2s, y2s, x3s, y3s, x4s, y4s);
    var d = basisToPoints(x1d, y1d, x2d, y2d, x3d, y3d, x4d, y4d);
    return multmm(d, adj(s));
}

function general2DProjectionFrom3points(
  x1s, y1s, x1d, y1d,
  x2s, y2s, x2d, y2d,
  x3s, y3s, x3d, y3d
) {
  var s = basisToPoints(x1s, y1s, x2s, y2s, x3s, y3s, x4s, y4s);
  var d = basisToPoints(x1d, y1d, x2d, y2d, x3d, y3d, x4d, y4d);
  return multmm(d, adj(s));
}

function project(m, x, y) {
    var v = multmv(m, [x, y, 1]);
    return [v[0] / v[2], v[1] / v[2]];
}

function solveAffineMatrix(r1, s1, t1, r2, s2, t2, r3, s3, t3) {
  const a = (((t2 - t3) * (s1 - s2)) - ((t1 - t2) * (s2 - s3))) / (((r2 - r3) * (s1 - s2)) - ((r1 - r2) * (s2 - s3)))
  const b = (((t2 - t3) * (r1 - r2)) - ((t1 - t2) * (r2 - r3))) / (((s2 - s3) * (r1 - r2)) - ((s1 - s2) * (r2 - r3)))
  const c = t1 - (r1 * a) - (s1 * b);
  return [a, b, c];
}

function deriveAffineTransform(a, b, c) {
  const e = 1e-15;
  a.xy.x -= e;
  a.xy.y += e;
  b.xy.x += e;
  b.xy.y -= e;
  c.xy.x += e;
  c.xy.y += e;
  x = solveAffineMatrix(
    a.xy.x, a.xy.y, a.latLonMeters.x,
    b.xy.x, b.xy.y, b.latLonMeters.x,
    c.xy.x, c.xy.y, c.latLonMeters.x,
  )
  y = solveAffineMatrix(
    a.xy.x, a.xy.y, a.latLonMeters.y,
    b.xy.x, b.xy.y, b.latLonMeters.y,
    c.xy.x, c.xy.y, c.latLonMeters.y,
  )
  return x.concat(y);
}

function CornersLatLonsFromThreePointsCoordsCal(width, height, gpsSeurantaCalString) {
  const calPtsRaw = gpsSeurantaCalString.split("|").map((x) => parseFloat(x));
  const proj = new SpheroidProjection();
  const calPts = [
    {
      latLonMeters: proj.latLonToMeters(new LatLon(calPtsRaw[1], calPtsRaw[0])),
      xy: new Point(calPtsRaw[2], calPtsRaw[3]),
    },
    {
      latLonMeters: proj.latLonToMeters(new LatLon(calPtsRaw[5], calPtsRaw[4])),
      xy: new Point(calPtsRaw[6], calPtsRaw[7]),
    },
    {
      latLonMeters: proj.latLonToMeters(new LatLon(calPtsRaw[9], calPtsRaw[8])),
      xy: new Point(calPtsRaw[10], calPtsRaw[11]),
    }
  ];
  const xyToLatLonMetersCoeffs = deriveAffineTransform(...calPts);
  function mapXYtoLatLon (xy) {
    const x = xy.x * xyToLatLonMetersCoeffs[0] + xy.y * xyToLatLonMetersCoeffs[1] + xyToLatLonMetersCoeffs[2];
    const y = xy.x * xyToLatLonMetersCoeffs[3] + xy.y * xyToLatLonMetersCoeffs[4] + xyToLatLonMetersCoeffs[5];
    return proj.metersToLatLon(new Point(x, y));
  }
  return [
    mapXYtoLatLon(new Point(0, 0)),
    mapXYtoLatLon(new Point(width, 0)),
    mapXYtoLatLon(new Point(width, height)),
    mapXYtoLatLon(new Point(0, height)),
  ];
}


function cornerCalTransform(width, height, top_left_latlon, top_right_latlon, bottom_right_latlon, bottom_left_latlon) {
    var proj = new SpheroidProjection();
    var top_left_meters = proj.latLonToMeters(top_left_latlon);
    var top_right_meters = proj.latLonToMeters(top_right_latlon);
    var bottom_right_meters = proj.latLonToMeters(bottom_right_latlon);
    var bottom_left_meters = proj.latLonToMeters(bottom_left_latlon);
    var matrix3d = general2DProjection(
        top_left_meters.x, top_left_meters.y, 0, 0,
        top_right_meters.x, top_right_meters.y, width, 0,
        bottom_right_meters.x, bottom_right_meters.y, width, height,
        bottom_left_meters.x, bottom_left_meters.y, 0, height
    )
    return function(latLon){
        var meters = proj.latLonToMeters(latLon);
        var xy = project(matrix3d, meters.x, meters.y);
        return new Point(xy[0], xy[1]);
    };
}

function cornerBackTransform(width, height, top_left_latlon, top_right_latlon, bottom_right_latlon, bottom_left_latlon) {
    var proj = new SpheroidProjection();
    var top_left_meters = proj.latLonToMeters(top_left_latlon);
    var top_right_meters = proj.latLonToMeters(top_right_latlon);
    var bottom_right_meters = proj.latLonToMeters(bottom_right_latlon);
    var bottom_left_meters = proj.latLonToMeters(bottom_left_latlon);
    var matrix3d = general2DProjection(
        0, 0, top_left_meters.x, top_left_meters.y,
        width, 0, top_right_meters.x, top_right_meters.y,
        width, height, bottom_right_meters.x, bottom_right_meters.y,
        0, height, bottom_left_meters.x, bottom_left_meters.y
    )
    return function(coords){
        var xy = project(matrix3d, coords.x, coords.y);
        return proj.metersToLatLon(new Point(xy[0], xy[1]));
    };
}

const dataURItoBlob = (dataURI) => {
  // convert base64/URLEncoded data component to raw binary data held in a string
  var byteString;
  if (dataURI.split(',')[0].indexOf('base64') >= 0)
      byteString = atob(dataURI.split(',')[1]);
  else
      byteString = unescape(dataURI.split(',')[1]);

  // separate out the mime component
  var mimeString = dataURI.split(',')[0].split(':')[1].split(';')[0];

  // write the bytes of the string to a typed array
  var ia = new Uint8Array(byteString.length);
  for (var i = 0; i < byteString.length; i++) {
      ia[i] = byteString.charCodeAt(i);
  }

  return new Blob([ia], {type:mimeString});
}

module.exports = {
    Point,
    LatLon,
    SpheroidProjection,
    cornerCalTransform,
    cornerBackTransform,
    dataURItoBlob,
    CornersLatLonsFromThreePointsCoordsCal
}