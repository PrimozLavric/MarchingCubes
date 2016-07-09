#include "Main.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <thread>
#include <fstream>
#include <chrono>
#include <math.h>

using namespace std;
using namespace std::chrono;

#define EXPORT_OBJ false

vec3<float> lerp(vec3<float> vec1, vec3<float> vec2, float alpha) {
	return vec1 + (vec2 - vec1) * alpha;
};

void marchingCubes(vector<vec3<float>> *vertices, const float *values, const vec3<int> &volDim, const int volZFull, const vec3<float> &voxDim, const float isoLevel, const int offset = 0) {

	// Actual position along edge weighted according to function values.
	vec3<float> vertList[12];


	// Calculate maximal possible axis value (used in vertice normalization)
	float maxX = voxDim.x * (volDim.x - 1);
	float maxY = voxDim.y * (volDim.y - 1);
	float maxZ = voxDim.z * (volZFull - 1);
	float maxAxisVal = max(maxX, max(maxY, maxZ));

	// Volume iteration
	for (int z = 0; z < volDim.z - 1; z++) {
		for (int y = 0; y < volDim.y - 1; y++) {
			for (int x = 0; x < volDim.x - 1; x++) {

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

				int p = x + (volDim.x * y) + (volDim.x * volDim.y * z),
					px = p + 1,
					py = p + volDim.x,
					pxy = py + 1,
					pz = p + volDim.x * volDim.y,
					pxz = px + volDim.x * volDim.y,
					pyz = py + volDim.x * volDim.y,
					pxyz = pxy + volDim.x * volDim.y;

				//							  X              Y                    Z
				vec3<float> position = { x * voxDim.x, y * voxDim.y, (z + offset) * voxDim.z };

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
				int bits = MC_EDGE_TABLE[cubeindex];

				// If no edge is triggered... skip
				if (bits == 0) continue;

				// Interpolate the positions based od voxel intensities
				float mu = 0.5;

				// bottom of the cube
				if (bits & 1) {
					mu = (isoLevel - value0) / (value1 - value0);
					vertList[0] = lerp(position, {position.x + voxDim.x, position.y, position.z}, mu);
				}
				if (bits & 2) {
					mu = (isoLevel - value1) / (value3 - value1);
					vertList[1] = lerp({position.x + voxDim.x, position.y, position.z}, {position.x + voxDim.x, position.y + voxDim.y, position.z}, mu);
				}
				if (bits & 4) {
					mu = (isoLevel - value2) / (value3 - value2);
					vertList[2] = lerp({ position.x, position.y + voxDim.y, position.z }, { position.x + voxDim.x, position.y + voxDim.y, position.z }, mu);
				}
				if (bits & 8) {
					mu = (isoLevel - value0) / (value2 - value0);
					vertList[3] = lerp(position, { position.x, position.y + voxDim.y, position.z }, mu);
				}
				// top of the cube
				if (bits & 16) {
					mu = (isoLevel - value4) / (value5 - value4);
					vertList[4] = lerp({ position.x, position.y, position.z + voxDim.z }, { position.x + voxDim.x, position.y, position.z + voxDim.z }, mu);
				}
				if (bits & 32) {
					mu = (isoLevel - value5) / (value7 - value5);
					vertList[5] = lerp({ position.x + voxDim.x, position.y, position.z + voxDim.z }, { position.x + voxDim.x, position.y + voxDim.y, position.z + voxDim.z }, mu);
				}
				if (bits & 64) {
					mu = (isoLevel - value6) / (value7 - value6);
					vertList[6] = lerp({ position.x, position.y + voxDim.y, position.z + voxDim.z }, { position.x + voxDim.x, position.y + voxDim.y, position.z + voxDim.z }, mu);
				}
				if (bits & 128) {
					mu = (isoLevel - value4) / (value6 - value4);
					vertList[7] = lerp({ position.x, position.y, position.z + voxDim.z }, { position.x, position.y + voxDim.y, position.z + voxDim.z }, mu);
				}
				// vertical lines of the cube
				if (bits & 256) {
					mu = (isoLevel - value0) / (value4 - value0);
					vertList[8] = lerp(position, { position.x, position.y, position.z + voxDim.z }, mu);
				}
				if (bits & 512) {
					mu = (isoLevel - value1) / (value5 - value1);
					vertList[9] = lerp({ position.x + voxDim.x, position.y, position.z }, { position.x + voxDim.x, position.y, position.z + voxDim.z }, mu);
				}
				if (bits & 1024) {
					mu = (isoLevel - value3) / (value7 - value3);
					vertList[10] = lerp({ position.x + voxDim.x, position.y + voxDim.y, position.z }, { position.x + voxDim.x, position.y + voxDim.y, position.z + voxDim.z }, mu);
				}
				if (bits & 2048) {
					mu = (isoLevel - value2) / (value6 - value2);
					vertList[11] = lerp({ position.x, position.y + voxDim.y, position.z }, { position.x, position.y + voxDim.y, position.z + voxDim.z }, mu);
				}

				// construct triangles -- get correct vertices from triTable.
				int i = 0;
				// "Re-purpose cubeindex into an offset into triTable."
				cubeindex <<= 4;

				while (MC_TRI_TABLE[cubeindex + i] != -1) {
					int index1 = MC_TRI_TABLE[cubeindex + i];
					int index2 = MC_TRI_TABLE[cubeindex + i + 1];
					int index3 = MC_TRI_TABLE[cubeindex + i + 2];

					// Add triangles vertices normalized with the maximal possible value
					(*vertices).push_back(vertList[index3] / maxAxisVal - 0.5);  // x
					(*vertices).push_back(vertList[index2] / maxAxisVal - 0.5);  // x
					(*vertices).push_back(vertList[index1] / maxAxisVal - 0.5);  // x

					i += 3;
				}
			}
		}
	}
};

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "Usage: script nThreadsMin nThreadsMax iterations\n";
		return 1;
	}

	vec3<int> size = { 512, 512, 512 };
	float axisMin = -10;
	float axisMax = 10;
	float axisRange = axisMax - axisMin;

	vector<float> scalarField;

	for (int k = 0; k < size.x; k++) {
		for (int j = 0; j < size.y; j++) {
			for (int i = 0; i < size.z; i++) {
				// actual values
				float x = axisMin + axisRange * i / (size.x - 1);
				float y = axisMin + axisRange * j / (size.y - 1);
				float z = axisMin + axisRange * k / (size.z - 1);
				scalarField.push_back(x*x + y*y - z*z - 25);
			
			}
		}
	}

	int nThreadsMin = atoi(argv[1]);
	int nThreadsMax = atoi(argv[2]);
	int iterations =  atoi(argv[3]);
	
	for (int nThreads = nThreadsMin; nThreads <= nThreadsMax; nThreads++) {

		vector<double> times;
		float sumTime = 0, avgTime;

		for (int it = 0; it < iterations; it++) {
			// TIMER
			auto start = system_clock::now();

			// Threads vector
			vector<thread> threads;
			// Results vector
			vector<vector<vec3<float>>> results(nThreads);

			// Thread work distribution
			int remainder = size.z % nThreads;
			int segment = size.z / nThreads;

			// First element pointer offset
			int offset = 0;

			// Z axis offset for vertice position calculation
			int zAxisOffset = 0;

			for (int i = 0; i < nThreads; i++) {
				// Distribute remainder among first (remainder) threads
				int segmentSize = (remainder-- > 0) ? segment + 1 : segment;

				// Padding needs to be added to correctly close the gaps between segments
				int paddedSegmentSize = (i != nThreads - 1) ? segmentSize + 1 : segmentSize;

				// Execute
				threads.push_back(thread(marchingCubes, &results[i], &scalarField[offset], vec3<int>{ size.x, size.y, paddedSegmentSize }, size.z, vec3<float>{ 1, 1, 1 }, 0.5, zAxisOffset));
				
				// Correct offsets for next iteration
				zAxisOffset += segmentSize;
				offset += segmentSize * size.x * size.y;
			}

			for (int i = 0; i < threads.size(); i++) {
				threads[i].join();
			}

			times.push_back(duration_cast<duration<double>>(system_clock::now() - start).count());
			sumTime += times.back();
		}

		// Average time
		avgTime = sumTime / iterations;

		// Standard deviation
		double sd = 0; 
		
		for (int i = 0; i < iterations; i++) {
			sd += pow((times[i] - avgTime), 2);
		}

		sd = sqrt(sd / iterations);
		

		cout << "Threads: " << nThreads << endl;
		cout << "Iterations: " << iterations << endl;
		cout << "Average time: " << avgTime << "s" << endl;
		cout << "Standard deviation: " << sd << endl;
		cout << "-------------------------------------------------------------" << endl;
	}
	system("pause");
	return 0;
}

