import java.util.ArrayList;

import java.lang.Math;

public class Main {

    static float[] lerp(float[] vec1, float[] vec2, float alpha){
        return new float[]{vec1[0] + (vec2[0] - vec1[0]) * alpha, vec1[1] + (vec2[1] - vec1[1]) * alpha, vec1[2] + (vec2[2] - vec1[2]) * alpha};
    }

    static void marchingCubes(float[] values, int[] volDim, int volZFull, float[] voxDim, float isoLevel, int offset, CallbackMC callback) {

        ArrayList<float[]> vertices = new ArrayList<>();
        // Actual position along edge weighted according to function values.
        float vertList[][] = new float[12][3];


        // Calculate maximal possible axis value (used in vertice normalization)
        float maxX = voxDim[0] * (volDim[0] - 1);
        float maxY = voxDim[1] * (volDim[1] - 1);
        float maxZ = voxDim[2] * (volZFull - 1);
        float maxAxisVal = Math.max(maxX, Math.max(maxY, maxZ));

        // Volume iteration
        for (int z = 0; z < volDim[2] - 1; z++) {
            for (int y = 0; y < volDim[1] - 1; y++) {
                for (int x = 0; x < volDim[0] - 1; x++) {

                    // Indices pointing to cube vertices
                    //              pyz  ___________________  pxyz
                    //                  /|                 /|
                    //                 / |                / |
                    //                /  |               /  |
                    //          pz   /___|______________/pxz|
                    //              |    |              |   |
                    //              |    |              |   |
                    //              | py |______________|___| pxy
                    //              |   /               |   /
                    //              |  /                |  /
                    //              | /                 | /
                    //              |/__________________|/
                    //             p                     px

                    int p = x + (volDim[0] * y) + (volDim[0] * volDim[1] * (z + offset)),
                            px = p + 1,
                            py = p + volDim[0],
                            pxy = py + 1,
                            pz = p + volDim[0] * volDim[1],
                            pxz = px + volDim[0] * volDim[1],
                            pyz = py + volDim[0] * volDim[1],
                            pxyz = pxy + volDim[0] * volDim[1];

                    //							  X              Y                    Z
                    float position[] = new float[]{x * voxDim[0], y * voxDim[1], (z + offset) * voxDim[2]};

                    // Voxel intensities
                    float value0 = values[p],
                            value1 = values[px],
                            value2 = values[py],
                            value3 = values[pxy],
                            value4 = values[pz],
                            value5 = values[pxz],
                            value6 = values[pyz],
                            value7 = values[pxyz];

                    // Voxel is active if its intensity is above isolevel
                    int cubeindex = 0;
                    if (value0 > isoLevel) cubeindex |= 1;
                    if (value1 > isoLevel) cubeindex |= 2;
                    if (value2 > isoLevel) cubeindex |= 8;
                    if (value3 > isoLevel) cubeindex |= 4;
                    if (value4 > isoLevel) cubeindex |= 16;
                    if (value5 > isoLevel) cubeindex |= 32;
                    if (value6 > isoLevel) cubeindex |= 128;
                    if (value7 > isoLevel) cubeindex |= 64;

                    // Fetch the triggered edges
                    int bits = TablesMC.MC_EDGE_TABLE[cubeindex];

                    // If no edge is triggered... skip
                    if (bits == 0) continue;

                    // Interpolate the positions based od voxel intensities
                    float mu = 0.5f;

                    // bottom of the cube
                    if ((bits & 1) != 0) {
                        mu = (isoLevel - value0) / (value1 - value0);
                        vertList[0] = lerp(position, new float[]{position[0] + voxDim[0], position[1], position[2]}, mu);
                    }
                    if ((bits & 2) != 0) {
                        mu = (isoLevel - value1) / (value3 - value1);
                        vertList[1] = lerp(new float[]{position[0] + voxDim[0], position[1], position[2]}, new float[]{position[0] + voxDim[0], position[1] + voxDim[1], position[2]}, mu);
                    }
                    if ((bits & 4) != 0) {
                        mu = (isoLevel - value2) / (value3 - value2);
                        vertList[2] = lerp(new float[]{position[0], position[1] + voxDim[1], position[2]}, new float[]{position[0] + voxDim[0], position[1] + voxDim[1], position[2]}, mu);
                    }
                    if ((bits & 8) != 0) {
                        mu = (isoLevel - value0) / (value2 - value0);
                        vertList[3] = lerp(position, new float[]{position[0], position[1] + voxDim[1], position[2]}, mu);
                    }
                    // top of the cube
                    if ((bits & 16) != 0) {
                        mu = (isoLevel - value4) / (value5 - value4);
                        vertList[4] = lerp(new float[]{position[0], position[1], position[2] + voxDim[2]}, new float[]{position[0] + voxDim[0], position[1], position[2] + voxDim[2]}, mu);
                    }
                    if ((bits & 32) != 0) {
                        mu = (isoLevel - value5) / (value7 - value5);
                        vertList[5] = lerp(new float[]{position[0] + voxDim[0], position[1], position[2] + voxDim[2]}, new float[]{position[0] + voxDim[0], position[1] + voxDim[1], position[2] + voxDim[2]}, mu);
                    }
                    if ((bits & 64) != 0) {
                        mu = (isoLevel - value6) / (value7 - value6);
                        vertList[6] = lerp(new float[]{position[0], position[1] + voxDim[1], position[2] + voxDim[2]}, new float[]{position[0] + voxDim[0], position[1] + voxDim[1], position[2] + voxDim[2]}, mu);
                    }
                    if ((bits & 128) != 0) {
                        mu = (isoLevel - value4) / (value6 - value4);
                        vertList[7] = lerp(new float[]{position[0], position[1], position[2] + voxDim[2]}, new float[]{position[0], position[1] + voxDim[1], position[2] + voxDim[2]}, mu);
                    }
                    // vertical lines of the cube
                    if ((bits & 256) != 0) {
                        mu = (isoLevel - value0) / (value4 - value0);
                        vertList[8] = lerp(position, new float[]{position[0], position[1], position[2] + voxDim[2]}, mu);
                    }
                    if ((bits & 512) != 0) {
                        mu = (isoLevel - value1) / (value5 - value1);
                        vertList[9] = lerp(new float[]{position[0] + voxDim[0], position[1], position[2]}, new float[]{position[0] + voxDim[0], position[1], position[2] + voxDim[2]}, mu);
                    }
                    if ((bits & 1024) != 0) {
                        mu = (isoLevel - value3) / (value7 - value3);
                        vertList[10] = lerp(new float[]{position[0] + voxDim[0], position[1] + voxDim[1], position[2]}, new float[]{position[0] + voxDim[0], position[1]+ voxDim[1], position[2] + voxDim[2]}, mu);
                    }
                    if ((bits & 2048) != 0) {
                        mu = (isoLevel - value2) / (value6 - value2);
                        vertList[11] = lerp(new float[]{position[0], position[1] + voxDim[1], position[2]}, new float[]{position[0], position[1] + voxDim[1], position[2] + voxDim[2]}, mu);
                    }

                    // construct triangles -- get correct vertices from triTable.
                    int i = 0;
                    // "Re-purpose cubeindex into an offset into triTable."
                    cubeindex <<= 4;

                    while (TablesMC.MC_TRI_TABLE[cubeindex + i] != -1) {
                        int index1 = TablesMC.MC_TRI_TABLE[cubeindex + i];
                        int index2 = TablesMC.MC_TRI_TABLE[cubeindex + i + 1];
                        int index3 = TablesMC.MC_TRI_TABLE[cubeindex + i + 2];

                        // Add triangles vertices normalized with the maximal possible value
                        vertices.add(new float[] {vertList[index3][0] / maxAxisVal - 0.5f, vertList[index3][1] / maxAxisVal - 0.5f, vertList[index3][2] / maxAxisVal - 0.5f});
                        vertices.add(new float[] {vertList[index2][0] / maxAxisVal - 0.5f, vertList[index2][1] / maxAxisVal - 0.5f, vertList[index2][2] / maxAxisVal - 0.5f});
                        vertices.add(new float[] {vertList[index1][0] / maxAxisVal - 0.5f, vertList[index1][1] / maxAxisVal - 0.5f, vertList[index1][2] / maxAxisVal - 0.5f});

                        i += 3;
                    }
                }
            }
        }

        callback.setVertices(vertices);
        callback.run();
    }

    public static void main(String[] args) {
        final int[] size = new int[]{256, 256, 256};

        final float[] scalarField = new float[size[0] * size[1] * size[2]];
        float axisMin = -10;
        float axisMax = 10;
        float axisRange = axisMax - axisMin;

        for (int k = 0; k < size[0]; k++) {
            for (int j = 0; j < size[1]; j++) {
                for (int i = 0; i < size[2]; i++) {
                    // actual values
                    float x = axisMin + axisRange * i / (size[0] - 1);
                    float y = axisMin + axisRange * j / (size[1] - 1);
                    float z = axisMin + axisRange * k / (size[2] - 1);
                    scalarField[k + size[1] * (j + size[2] * i)] = (x * x + y * y - z * z - 25);
                }
            }
        }

        int nThreadsMin = Integer.parseInt(args[0]);
        int nThreadsMax = Integer.parseInt(args[1]);
        int iterations = Integer.parseInt(args[2]);


        for (int nThreads = nThreadsMin; nThreads <= nThreadsMax; nThreads++) {

            final ArrayList<Double> times = new ArrayList<>();

            for (int it = 0; it < iterations; it++) {

                // TIMER
                final long start = System.currentTimeMillis();

                ArrayList<Thread> threads = new ArrayList<>();
                final ArrayList<ArrayList<float []>> results = new ArrayList<>();


                // Thread work distribution
                int remainder = size[2] % nThreads;
                int segment = size[2] / nThreads;

                // Z axis offset for vertice position calculation
                int zAxisOffset = 0;


                for (int i = 0; i < nThreads; i++) {
                    // Distribute remainder among first (remainder) threads
                    int segmentSize = (remainder-- > 0) ? segment + 1 : segment;

                    // Padding needs to be added to correctly close the gaps between segments
                    final int paddedSegmentSize = (i != nThreads - 1) ? segmentSize + 1 : segmentSize;


                    // Finished callback
                    final CallbackMC callback = new CallbackMC() {
                        @Override
                        public void run() {
                            results.add(getVertices());
                        }
                    };

                    // Java...
                    final int finalZAxisOffset = zAxisOffset;

                    // Start the thread
                    Thread t = new Thread() {
                        public void run() {
                            marchingCubes(scalarField, new int[]{size[0], size[1], paddedSegmentSize}, size[2], new float[]{1.0f, 1.0f, 1.0f}, 0.5f, finalZAxisOffset, callback);
                        }
                    };

                    threads.add(t);
                    t.start();

                    // Correct offsets for next iteration
                    zAxisOffset += segmentSize;
                }

                // Join the threads
                for (int i = 0; i <threads.size(); i ++) {
                    try {
                        threads.get(i).join();
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }

                // Time measurement
                long end= System.currentTimeMillis();
                times.add((end - start) / 1000.0);
            }

            double sumTime = 0.0, avgTime;

            for (int i = 0; i < times.size(); i++) {
                sumTime += times.get(i);
            }
            // Average time
            avgTime = sumTime / iterations;

            // Standard deviation
            double sd = 0;

            for (int i = 0; i < times.size(); i++) {
                sd += Math.pow((times.get(i) - avgTime), 2);
            }

            sd = Math.sqrt(sd / iterations);


            System.out.println("Threads: " + nThreads);
            System.out.println("Iterations: " + iterations);
            System.out.println("Average time: " + avgTime + "s");
            System.out.println("Standard deviation: " + sd);
            System.out.println("-------------------------------------------------------------");
        }


    }
}
