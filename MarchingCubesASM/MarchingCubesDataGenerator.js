/**
 * Created by Primoz on 7.6.2016.
 */

onmessage = function(msg) {
    // number of cubes along a side
    var size = msg.data.size;
    var planeEq = msg.data.planeEq;

    var axisMax = 10, axisMin = - 10;
    var axisRange = (axisMax - axisMin);

    var values = new Float32Array(size.x * size.y * size.z);

    var idx = 0;
    // Generate a list of 3D positions and values
    for (var k = 0; k < size.z; k++) {
        for (var j = 0; j < size.y; j++) {
            for (var i = 0; i < size.x; i++) {
                // actual values
                var x = axisMin + axisRange * i / (size.x - 1);
                var y = axisMin + axisRange * j / (size.y - 1);
                var z = axisMin + axisRange * k / (size.z - 1);

                var value = eval(planeEq);
                values[idx] = value;

                idx ++;
            }
        }
    }

    postMessage(values.buffer, [values.buffer]);
};