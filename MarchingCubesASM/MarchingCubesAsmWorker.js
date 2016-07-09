/**
 * Created by Primoz on 6.6.2016.
 */
self.importScripts('asmMC.js');


var msgCount = 0;

var dimensions;
var voxelDim;
var isoLevel;
var offset;

var values;

var jsBridge = null;

onmessage = function(msg) {
    // First message should specify the meta data
    if (msgCount <= 0) {
        var meta = msg.data;
        dimensions = meta.dimensions;
        voxelDim = meta.voxelDimensions;
        isoLevel = meta.isoLevel;
        // If parallelized take correctly offset z axis
        offset = (meta.dimensions.offset) ? meta.dimensions.offset : 0;

        msgCount++;
        return;
    }
    else if (msgCount === 1) {
        values = new Float32Array(msg.data);
    }

    if(!jsBridge)
        jsBridge = new JsBridge();


    var vertices = jsBridge.marchingCubes(dimensions.x, dimensions.y, dimensions.z, dimensions.zFull, voxelDim.x, voxelDim.y, voxelDim.z, isoLevel, offset, values);
    
    postMessage(vertices.buffer, [vertices.buffer]);
};