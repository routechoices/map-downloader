#! /usr/bin/env node 
const { loadImage } = require('@napi-rs/canvas');
const url = require('url')
const fetch = require('node-fetch')
const express = require('express')
const stream = require('stream')
const { drawRoute, saveKMZ, drawMapWithCourse, getProj4Def } = require('./helpers')
const sharp = require('sharp');
const fileUpload = require('express-fileupload');
const { readOcad } = require('ocad2geojson')
const OcadTiler = require('ocad-tiler')
const proj4 = require('proj4')
const { render } = require('./ocad_render')

const app = express()

const getLiveloxMap = async (req, res, next) => {
    const liveloxUrl = req.body.url
    if (!liveloxUrl.startsWith('https://www.livelox.com/')) { 
        return res.status(400).send('invalid url domain')
    }
    let classId = ''
    try {
        classId = url.parse(liveloxUrl, true).query.classId
    } catch (e) {
        return res.status(400).send('no class id provided' )
    }
    let data = {}
    try {
        const res = await fetch("https://www.livelox.com/Data/ClassInfo", {
            "headers": {
                "accept": "application/json",
                "content-type": "application/json",
                "X-Requested-With": "XMLHttpRequest",
            },
            "body": JSON.stringify({
                "classIds":[classId],
                "courseIds":null,
                "relayLegs":[],
                "relayLegGroupIds":[],
                "includeMap":true,
                "includeCourses":true,
                "skipStoreInCache":false
            }),
            "method": "POST"
        });
        data = await res.json()
    } catch (e) {
        return res.status(400).send('could not reach livelox server')
    }
    const eventData = data.general
    const blobUrl = eventData?.classBlobUrl
    if (!blobUrl) {
        return res.status(400).send('cannot not figure blob url')
    }
    let blobData = null
    try {
        const res = await fetch(blobUrl, {
            "headers": {
                "accept": "application/json",
                "content-type": "application/json",
                "X-Requested-With": "XMLHttpRequest",
            }
        });
        blobData = await res.json()
    } catch (e) {
        return res.status(400).send('could not reach blob url')
    }

    let mapUrl, mapBound, mapResolution, route, mapName
    try {
        mapData = blobData.map
        mapUrl = mapData.url;
        mapBound = mapData.boundingQuadrilateral.vertices
        mapResolution = mapData.resolution
        route = blobData.courses.map((c) => c.controls)
        mapName = mapData.name
    } catch (e) {
        return res.status(400).send('could not parse livelox data')
    }
    try {
        const mapImg = await loadImage(mapUrl)
        const [outCanvas, bounds] = drawRoute(mapImg, mapBound, route, mapResolution)
        const imgBlob = outCanvas.toBuffer('image/png')
        const outImgBlob = await sharp(imgBlob).webp().toBuffer()
        let buffer
        let mime
        let filename
        if (!req.body.type || req.body.type === 'webp') {
            buffer = outImgBlob
            mime = 'image/webp'
            filename = `${mapName}_${bounds[3].lat}_${bounds[3].lon}_${bounds[2].lat}_${bounds[2].lon}_${bounds[1].lat}_${bounds[1].lon}_${bounds[0].lat}_${bounds[0].lon}_.webp`
        } else if(req.body.type === 'kmz') {
            buffer = await saveKMZ(
                mapName,
                {
                    top_left: bounds[3],
                    top_right: bounds[2],
                    bottom_right: bounds[1],
                    bottom_left: bounds[0]
                },
                outImgBlob
            )
            mime = 'application/vnd.google-earth.kmz'
            filename = `${mapName}.kmz`
        } else {
            return res.status(400).send('invalid type' )
        }
        var readStream = new stream.PassThrough()
        readStream.end(buffer)
        res.set('Content-disposition', 'attachment; filename="' + filename.replace(/\\/g, '_').replace(/"/g, '\\"') + '"')
        res.set('Content-Type', mime)
        readStream.pipe(res)
    } catch (e) {
        return res.status(500).send('Something went wrong... '+e.message)
    }
}

const getRGClasses = async (req, res, next) => {
    const eventUrl = req.body.url;
    const parsedUrl = url.parse(eventUrl, true)
    const gadgetRootPath = parsedUrl.pathname.split('/').slice(0, -2).join('/') + "/kartat"
    const eventId = parseInt(parsedUrl.query.id, 10);
    const dataFile = parsedUrl.protocol + '//' + parsedUrl.host + '/' + gadgetRootPath + "/sarjat_" + eventId + ".txt";
    console.log(dataFile)
    const classesFileRequest = await fetch(dataFile)
    if (classesFileRequest.status != 200) {
      return res.status(200).send({error: "Cannot access classes file"})
    }
    const classesFile = await classesFileRequest.text();
    console.log(classesFile)
    const lines = classesFile.split('\n').map((l) => l.trim()).filter(Boolean)
    const classes = lines.map((line) => {
      const data = line.split('|');
      return [data[0], data.slice(1).join('|')]
    }).filter(Boolean)
    return res.status(200).send({classes, eventUrl: eventUrl})
  }

const getRGMap = async (req, res, next) => {
    const eventUrl = req.body.url;
    console.log(eventUrl)
    if (!eventUrl  || !req.body.classId) {
      return res.status(200).send({error: "Missing parameters"})
    }
    const parsedUrl = url.parse(eventUrl, true)
    const gadgetRootPath = parsedUrl.pathname.split('/').slice(0, -2).join('/') + "/kartat"
    const eventId = parseInt(parsedUrl.query.id, 10);
  
    const cFileUrl = parsedUrl.protocol + '//' + parsedUrl.host + '/' + gadgetRootPath + "/kilpailijat_" + eventId + ".txt";
    console.log(cFileUrl)
    const competitorFileRequest = await fetch(cFileUrl)
    if (competitorFileRequest.status != 200) {
      return res.status(200).send({error: "Cannot access competitor file"})
    }
    const cFile = await competitorFileRequest.text();
    const clines = cFile.split('\n').map((l) => l.trim()).filter(Boolean)
    const routesIdsRaw = clines.map((line) => {
      return line.split('|');
    }).filter((d) => {
      return d?.[1] == req.body.classId
    }).map((d) => d?.[6] || d?.[1])
    const routesIds = [...new Set(routesIdsRaw)];
    if (!routesIds.length) {
      return res.status(200).send({error: "Cannot find routes in competitors file"})
    }
    const dataFile = parsedUrl.protocol + '//' + parsedUrl.host + '/' + gadgetRootPath + "/ratapisteet_" + eventId + ".txt";
    console.log(dataFile)
    const routesFileRequest = await fetch(dataFile)
    if (routesFileRequest.status != 200) {
      return res.status(200).send({error: "Cannot access routes file"})
    }
    const routesFile = await routesFileRequest.text();
    console.log(routesFile)
    const lines = routesFile.split('\n').map((l) => l.trim()).filter(Boolean)
    const routesData = lines.map((line) => {
      const data = line.split('|');
      return [data[0], data.slice(1).join('|')]
    }).filter((d) => {
      return routesIds.includes(d?.[0])
    }).map((d) => d?.[1])
    if (!routesData.length) {
      return res.status(200).send({error: "Cannot find routes in routes file"})
    }
    const coordinates = routesData.map(routeData => routeData.split('N').map((xy) => xy && xy.split(';').map((x) => parseInt(x, 10))).filter(Boolean))
    console.log(coordinates)
    
    const mapListFileURL = parsedUrl.protocol + '//' + parsedUrl.host + '/' + gadgetRootPath + "/kisat.txt";
    const mapFileRequest = await fetch(mapListFileURL)
  
    if (mapFileRequest.status != 200) {
      return res.status(200).send({error: "Cannot access routes file"})
    }
    const mapListFile = await mapFileRequest.text();
    console.log(mapListFile)
    const klines = mapListFile.split('\n').map((l) => l.trim()).filter(Boolean)
    const mapId = klines.map((line) => {
      return line.split('|');
    }).find((d) => {
      return d?.[0] == eventId
    })?.[1]
    console.log(mapId)
    
    const mapURL = parsedUrl.protocol + '//' + parsedUrl.host + '/' + gadgetRootPath + "/" + parseInt(mapId, 10) + ".jpg";
    const mapImg = await loadImage(mapURL)
    const resultImg = drawMapWithCourse(mapImg, coordinates)
    const buffer = resultImg.toBuffer('image/jpeg')
    const mime = 'image/jpeg'
    const filename = `map.jpg`
    var readStream = new stream.PassThrough()
    readStream.end(buffer)
    res.set('Content-disposition', 'attachment; filename="' + filename.replace(/\\/g, '_').replace(/"/g, '\\"') + '"')
    res.set('Content-Type', mime)
    readStream.pipe(res)
  }

  const getOcadMap = async (req, res, next) => {
    if(!req.files?.ocad_file) {
       return res.status(400).send('no file sent')
    }
    if (!['webp', 'kmz'].includes(req.body.type)){
       return res.status(400).send('invalid output format')
    }
    const uploadedFile = req.files.ocad_file;

    const ocadFile = await readOcad(uploadedFile.data);
    const mapCrs = ocadFile.getCrs();
    const mapGeoRef = mapCrs.code !== 0;
    if (!mapGeoRef) {
        return res.status(400).send("Map not geo-referenced")
    }

    let northEast, northWest, southWest, southEast;
    if (mapGeoRef) {
        const proj4Def = await getProj4Def(mapCrs.code);
        const projectedBounds = ocadFile.getBounds(mapCrs.toProjectedCoord.bind(mapCrs))

        proj4.defs('WGS84', "+title=WGS 84 (long/lat) +proj=longlat +ellps=WGS84 +datum=WGS84 +units=degrees");
        proj4.defs('OcadFile', proj4Def);

        northEast = proj4('OcadFile', 'WGS84', [projectedBounds[2], projectedBounds[3]])
        northWest = proj4('OcadFile', 'WGS84', [projectedBounds[0], projectedBounds[3]])
        southWest = proj4('OcadFile', 'WGS84', [projectedBounds[0], projectedBounds[1]])
        southEast = proj4('OcadFile', 'WGS84', [projectedBounds[2], projectedBounds[1]])
    }
    const tiler = new OcadTiler(ocadFile)
    tileBounds = tiler.bounds
    const imgBlob = await render(tiler, tileBounds, 1, {
        format: "webp",
        exportHidden: false,
        applyGrivation: false,
        fill: "#fff"
    })
    let filename = ""
    let mime = ""
    let out = null;

    if (!req.body.type || req.body.type === 'webp') {
        out = imgBlob
        mime = 'image/webp'
        if (mapGeoRef) {
            filename = uploadedFile.name.slice(0, -4) + "_" + northWest[1] + "_" + northWest[0] + "_" + northEast[1] + "_" + northEast[0] + "_" + southEast[1] + "_" + southEast[0] + "_" + southWest[1] + "_" + southWest[0] + "_" +".webp"
        } else {
            filename = uploadedFile.name.slice(0, -4) + ".webp"
        }
    } else if(req.body.type === 'kmz') {
        const mapName = uploadedFile.name.slice(0, -4)
        out = await saveKMZ(
            mapName,
            {
                top_left: {lat: northWest[1], lon: northWest[0]},
                top_right: {lat: northEast[1], lon: northEast[0]},
                bottom_right: {lat: southEast[1], lon: southEast[0]},
                bottom_left: {lat: southWest[1], lon: southWest[0]},
            },
            imgBlob
        )
        mime = 'application/vnd.google-earth.kmz'
        filename = `${mapName}.kmz`
    } else {
        return res.status(400).send('invalid type' )
    }

    const readStream = new stream.PassThrough()
    readStream.end(out)

    res.set('Content-Type', mime + '; charset=utf-8')
    res.set('Content-Disposition', "attachment; filename*=UTF-8''" + escape(filename))
    readStream.pipe(res)
}

app.use(express.urlencoded({extended: true}))
app.use(express.json())
app.use(fileUpload({
    createParentPath: true
}));
app.use(express.static('public'))

app.post('/api/get-livelox-map', getLiveloxMap)
app.post('/api/get-routegadget-classes', getRGClasses)
app.post('/api/get-routegadget-map', getRGMap)
app.post('/api/get-ocad-map', getOcadMap)

module.exports = app
