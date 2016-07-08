/**
 * Created by Primoz on 31.5.2016.
 */

MarchingCubes = class {

    constructor () {
        this._jobQueue = [];
        this._isRunning = false;
    }

    extractMesh (meta, values, nThreads, onLoad) {
        if (!onLoad) {
            console.error("Tried to execute marching cubes without onLoad callback!");
            return;
        }

        // Add the job in the jobQueue
        this._jobQueue.push({meta: meta, values: values, nThreads: nThreads, onLoad: onLoad});

        // If no job is currently executing. Start the execution.
        if (!this._isRunning) {
            this._executeNextJob();
        }
    }

    _executeNextJob() {
        this._jobQueue[0].values = deepCopy(this._jobQueue[0].values);

        var start = new Date(), stop;
        // Set the job running flag
        this._isRunning = true;
        var self = this;

        // Split the work among workers
        var meta = this._jobQueue[0].meta;
        var nThreads = this._jobQueue[0].nThreads;

        // Array for combined results finish counters
        var vertices = [];

        // SINGLE WORKER EXECUTION
        if (nThreads <= 1 || meta.dimensions.z < nThreads) {
            let worker = new Worker("./MarchingCubesWorker.js");


            // When single worker is used.. When the result message comes.. immediately execute the callback and move to the next task
            worker.onmessage = function(result) {
                if (result.data instanceof ArrayBuffer) {
                    // Data received
                    vertices.push(new Float32Array(result.data));
                }
                else if (result.data.type === "finished") {

                    // Last (finish) message contains length of the last buffer
                    var lastLength = result.data.lastLength;
                    var lastBuffer = vertices[vertices.length - 1];
                    var clampedBuffer = new Float32Array(lastLength);

                    // Copy element and replace last buffer with clamped buffer
                    clampedBuffer.set(lastBuffer.slice(0, lastLength));
                    vertices[vertices.length - 1] = clampedBuffer;
                    lastBuffer = null;

                    // Terminate worker once finished
                    worker.terminate();

                    // Performance measurement
                    stop = new Date();
                    var time = (stop - start) / 1000;

                    // Notify user about the results
                    self._jobQueue[0].onLoad(vertices, time);
                    self._jobQueue.shift();

                    // Check if there is anything left in the queue
                    if (self._jobQueue.length !== 0) {
                        self._executeNextJob();
                    }
                    else {
                        self._isRunning = false;
                    }
                }
            };


            // Start the worker task
            // Pass meta data
            worker.postMessage({dimensions: meta.dimensions, voxelDimensions: meta.voxelDimensions, isoLevel: meta.isoLevel, valuesType: this._jobQueue[0].values.constructor.name});

            // Pass data
            worker.postMessage(this._jobQueue[0].values.buffer, [this._jobQueue[0].values.buffer]);
        }
        // MULTIPLE WORKERS EXECUTION
        else {
            // Calculate segment sizes (work distribution)
            var remainder = meta.dimensions.z % nThreads;
            var segment = Math.trunc(meta.dimensions.z / nThreads);

            // Counts number of threads that finished
            var finishedCounter = 0;

            // Work segmentation offsets
            var offset = 0;
            var zAxisOffset = 0;

            // Error management
            var workers = [];

            for (let i = 0; i < nThreads; i++) {
                // Correctly distribute the remainder
                let size = (remainder-- > 0) ? segment + 1 : segment;
                // Padding needs to be added to correctly close the gaps between segments
                let paddedSize = (i !== nThreads - 1) ? size + 1 : size;
                let chunkSize = paddedSize * meta.dimensions.x * meta.dimensions.y;

                // Split the data (slice makes shallow copy)
                let valuesSegment = this._jobQueue[0].values.slice(offset, offset + chunkSize);

                offset += size * meta.dimensions.x * meta.dimensions.y;

                // Initialize and start workers
                let worker = new Worker("./MarchingCubesWorker.js");
                workers.push(worker);

                worker.onmessage = function (result) {
                    if (result.data instanceof ArrayBuffer) {
                        // New vertice batch
                        vertices.push(new Float32Array(result.data));
                    }
                    else if (result.data.type === "finished") {
                        // Worker finished
                        // Last (finish) message contains length of the last buffer
                        var lastLength = result.data.lastLength;
                        var lastBuffer = vertices[vertices.length - 1];
                        var clampedBuffer = new Float32Array(lastLength);

                        // Copy element and replace last buffer with clamped buffer
                        clampedBuffer.set(lastBuffer.slice(0, lastLength));
                        vertices[vertices.length - 1] = clampedBuffer;
                        lastBuffer = null;

                        // Thread finished counter
                        finishedCounter ++;

                        // Clean up
                        worker.terminate();

                        // When the last worker finishes.. return the combined result via callback
                        if (finishedCounter === nThreads) {

                            stop = new Date();
                            var time = (stop - start) / 1000;

                            // Notify user about the results
                            self._jobQueue[0].onLoad(vertices, time);
                            self._jobQueue.shift();

                            // Check if there is anything left in the queue
                            if (self._jobQueue.length !== 0) {
                                self._executeNextJob();
                            }
                            else {
                                self._isRunning = false;
                            }
                        }
                    }
                };


                // Pass meta data
                worker.postMessage({dimensions: {x: meta.dimensions.x, y: meta.dimensions.y, z: paddedSize, zFull: meta.dimensions.z, offset: zAxisOffset},
                    voxelDimensions: meta.voxelDimensions, isoLevel: meta.isoLevel, valuesType: this._jobQueue[0].values.constructor.name});
                // Pass data
                worker.postMessage(valuesSegment.buffer, [valuesSegment.buffer]);

                // Update z axis offset
                zAxisOffset += size;
            }
        }
        this._jobQueue[0].values = null;
    }
};

function deepCopy(array) {
    var deepCopy = eval("new " + array.constructor.name + "(array.length);");

    for (var i = 0; i < array.length; i++) {
        deepCopy[i] = array[i];
    }

    return deepCopy;
}

