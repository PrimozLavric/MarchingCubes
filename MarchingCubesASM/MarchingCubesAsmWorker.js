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

// TODO: Implement support for other value types

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
        values = new Uint8Array(msg.data);
    }

    // VALUES MEMORY ALLOCATION
    // Get positions byte size and alloc space on the heap
    var nValBytes = values.length * values.BYTES_PER_ELEMENT;
    var valPtr = Module._malloc(nValBytes);

    // Copy values data to the heap
    var valHeapAlloc = new Uint8Array(Module.HEAPU8.buffer, valPtr, nValBytes);
    valHeapAlloc.set(values);


    // Returns pointer to result
    var rezPoint = Module.ccall('marchingCubes',
        'number', ['number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number'],
        [dimensions.x, dimensions.y, dimensions.z, dimensions.zFull, voxelDim.x, voxelDim.y, voxelDim.z, isoLevel, offset, valHeapAlloc.byteOffset]);

    Module._free(valPtr);

    // Parse results
    var rezSize = getValue(rezPoint, 'float');

    console.log(rezSize);
    var vertices = new Float32Array(rezSize-1);
    for (var i = 0; i < rezSize; i++) {
        vertices[i] = getValue(rezPoint + (4 * (i+1)), 'float')
    }

    console.log(vertices.length);

    postMessage(vertices.buffer, [vertices.buffer]);
};