#include "Main.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <thread>
#include <fstream>
#include <chrono>
#include <math.h>
#include <string>
#include <regex>

using namespace std;
using namespace std::chrono;

vec3<float> lerp(vec3<float> vec1, vec3<float> vec2, float alpha) {
	return vec1 + (vec2 - vec1) * alpha;
};

template<typename T>
void marchingCubes(vector<vec3<float>> *vertices, const T *values, const vec3<int> &volDim, const int volZFull, const vec3<float> &voxDim, const T isoLevel, const int offset = 0) {

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
				T	value0 = values[p],
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

// Generates volume with hyperboloid equation: x^2 + y^2 - z^2 - 25 with axes range [-10, 10].
template<typename T>
vector<T> generateScalarVolume(const vec3<int> &size) {
	float axisMin = -10;
	float axisMax = 10;
	float axisRange = axisMax - axisMin;

	vector<T> scalarField;

	for (int k = 0; k < size.x; k++) {
		for (int j = 0; j < size.y; j++) {
			for (int i = 0; i < size.z; i++) {
				// actual values
				T x = axisMin + axisRange * i / (size.x - 1);
				T y = axisMin + axisRange * j / (size.y - 1);
				T z = axisMin + axisRange * k / (size.z - 1);
				scalarField.push_back(x*x + y*y - z*z - 25);

			}
		}
	}

	return scalarField;
}
 
template<typename T>
void extractInputHandler(ifstream &inputFileStream, ofstream &outFileStream, const vec3<int> &size, const vec3<float> &voxSize, const T &isoValue, const int &nThreads) {
	T *scalarField;
	vector<T> helperVec;

	// Check if input file was specified
	if (inputFileStream.is_open()) {
		cout << "PROGRESS: Reading input data." << endl;
		// Read the data from the input file
		scalarField = new T[size.x * size.z * size.y];
		long fileSize = size.x * size.z * size.y * sizeof(T);
		
		inputFileStream.read((char*)scalarField, fileSize);
		inputFileStream.close();
	}
	else {
		cout << "PROGRESS: Generating volume data." << endl;
		helperVec = generateScalarVolume<T>(size);
		scalarField = &helperVec[0];
	}

	cout << "PROGRESS: Executing marching cubes." << endl;

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
		threads.push_back(thread(marchingCubes<T>, &results[i], scalarField + offset, vec3<int>{ size.x, size.y, paddedSegmentSize }, size.z, voxSize, isoValue, zAxisOffset));

		// Correct offsets for next iteration
		zAxisOffset += segmentSize;
		offset += segmentSize * size.x * size.y;
	}

	for (int i = 0; i < threads.size(); i++) {
		threads[i].join();
	}

	cout << "PROGRESS: Writing results to file." << endl;

	// WRITE .OBJ TO FILE
	if (!outFileStream.is_open()) {
		cout << "Error: File output stream is closed!";
		return;
	}

	int idx = 0;
	for (int i = 0; i < results.size(); i++) {
		for (int j = 0; j < results[i].size(); j++) {
			outFileStream << "v " << results[i][j].x << " " << results[i][j].y << " " << results[i][j].z << endl;
			
			if (idx % 3 == 0) {
				outFileStream << "f " << idx + 1 << " " << idx + 2 << " " << idx + 3 << endl;
				outFileStream.flush();
			}

			idx += 1;
		}
	}
		
	outFileStream.close();
}

template<typename T>
void benchmarkHandler(ifstream &inputFileStream, ofstream &outFileStream, const vec3<int> &size, const vec3<float> &voxSize, const T &isoValue, const int &nThreadsMin, const int &nThreadsMax, const int &iterations) {
	
	T *scalarField;
	vector<T> helperVec;

	// Check if input file was specified
	if (inputFileStream.is_open()) {
		cout << "PROGRESS: Reading input data." << endl;
		// Read the data from the input file
		scalarField = new T[size.x * size.z * size.y];
		long fileSize = size.x * size.z * size.y * sizeof(T);

		inputFileStream.read((char*)scalarField, fileSize);
		inputFileStream.close();
	}
	else {
		cout << "PROGRESS: Generating volume data." << endl;
		helperVec = generateScalarVolume<T>(size);
		scalarField = &helperVec[0];
	}

	cout << "PROGRESS: Starting benchmarks." << endl;
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
				threads.push_back(thread(marchingCubes<T>, &results[i], scalarField + offset, vec3<int>{ size.x, size.y, paddedSegmentSize }, size.z, voxSize, isoValue, zAxisOffset));

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

		// If otuput file is specified.. write the results to file
		if (outFileStream.is_open()) {
			outFileStream << "Threads: " << nThreads << endl;
			outFileStream << "Iterations: " << iterations << endl;
			outFileStream << "Average time: " << avgTime << "s" << endl;
			outFileStream << "Standard deviation: " << sd << endl;
			outFileStream << "-------------------------------------------------------------" << endl;
		}
	}
	cout << "PROGRESS: Finished benchmarks." << endl;
}


int main(int argc, char* argv[]) {
	static string usage = "This script may be executed in either benchmark or extract mode. Mode is specified by the first parameter [benchmark, extract].\nParameters: \n\t-input-vol\t Specifies path to the input volume. If this parameter is set volume dimensions(-vol-dim), data type(-data-type) and iso value(-iso) must also be given.\n\t-vol-dim\t Specifies the generated/read volume dimensions. Dimensions should be given as unsigned integers in format; -vol-dim X Y Z.\n\t-data-type\t Specifies the input file or generated data type. Options [char, uchar, short, ushort, int, uint, float, double].\n\t-vox-dim\t Specifies voxel dimensions used in mesh construction. Dimensions should be given as floating point numbers in format: -vox-dim X Y Z.\n\t-nThread\t Number of threads used in Marching cubes algorithm.This parameter can be either given as a single unsigned integer value or two unsigned integer values in benchmark mode, specifying the range of thread executions that will be tested.\n\t-iter\t\t Used only in benchmark mode to determine how many iterations should be executed for each configuration.\n\t-iso\t\t Isovalue that is used as a threshold for determining active voxels. Type should match the data type.\n\t-o\t\t Path to output file. In extract mode the mesh is written to file in .obj format [required]. In benchmark mode the results are written to file.\n";
	static regex anyNumRegex("(([1-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?");
	static regex uintNumRegex("[1-9][0-9]*");

	if (argc < 2) {
		cout << usage << endl;
		return 1;
	}
	else if (strcmp(argv[1], "-help") == 0) {
		cout << usage << endl;
		return 0;
	}

	bool benchmark = false;

	int nThreadsMin = thread::hardware_concurrency();
	if (nThreadsMin == 0) {
		nThreadsMin = 1;
	}
	int nThreadsMax = nThreadsMin;

	ifstream inputFileStream;
	ofstream outFileStream;
	string type;
	string isoValueStr;
	int iterations = 10;	// Default 10 iterations per benchmark

	bool customSizeSpecified = false;
	vec3<int> size = { 64, 64, 64 };
	vec3<float> voxSize = { 1.0f, 1.0f, 1.0f };

	// START PARAMETER PARSING
	// Read execution type
	if (strcmp(argv[1], "benchmark") == 0) {
		benchmark = true;
	}
	else if (strcmp(argv[1], "extract") != 0) {
		cout << "Invalid execution type. Valid options [extract, benchmark]" << endl;
		return 1;
	}

	// Flag parsing
	for (int i = 2; i < argc; i++) {
		if (strcmp(argv[i], "-input-vol") == 0) {
			// Volume path specified
			// Output file path is specified
			if (i + 1 >= argc || argv[i + 1][0] == '-') {
				cout << "Missing file path after -input-vol flag." << endl;
				return 1;
			}

			// Store the file name and offset iterator
			inputFileStream = ifstream(argv[++i], ios::in | ios::binary);

			if (!inputFileStream.good()) {
				cout << "Specified volume file does not exist." << endl;
				return 1;
			}
		}
		else if (strcmp(argv[i], "-vol-dim") == 0) {
			// Volume dimensions are given
			if (i + 3 >= argc || argv[i + 1][0] == '-' || argv[i + 2][0] == '-' || argv[i + 3][0] == '-') {
				cout << "Missing volume dimensions after -vol-dim flag." << endl;
				return 1;
			}

			string x(argv[++i]);
			string y(argv[++i]);
			string z(argv[++i]);

			if (!regex_match(x, uintNumRegex) || !regex_match(y, uintNumRegex) || !regex_match(z, uintNumRegex)) {
				cout << "Invalid volume dimensions format. Specify dimensions as three unsigned integers." << endl;
				return 1;
			}

			customSizeSpecified = true;
			size = { atoi(x.c_str()), atoi(y.c_str()), atoi(z.c_str()) };
		}
		else if (strcmp(argv[i], "-vox-dim") == 0) {
			// Voxel dimensions are given
			if (i + 3 >= argc) {
				cout << "Missing voxel dimensions after -vox-dim flag." << endl;
				return 1;
			}

			string x(argv[++i]);
			string y(argv[++i]);
			string z(argv[++i]);

			if (!regex_match(x, anyNumRegex) || !regex_match(y, anyNumRegex) || !regex_match(z, anyNumRegex)) {
				cout << "Invalid voxel dimensions format. Specify voxel dimensions as three positive floats." << endl;
				return 1;
			}
 
			
			voxSize = { (float) atof(x.c_str()), (float) atof(y.c_str()), (float) atof(z.c_str()) };
		}
		else if (strcmp(argv[i], "-nThread") == 0) {
			// Number of threads is given
			// FIRST VALUE
			if (i + 1 >= argc || argv[i + 1][0] == '-') {
				cout << "Missing number or range of threads after -nThread flag." << endl;
				return 1;
			}

			// Validate first number
			string tmp(argv[++i]);

			if (!regex_match(tmp, uintNumRegex)) {
				cout << "Invalid nThread value format. Specify unsigned integer value or two if range." << endl;
				return 1;
			}

			// Parse C-str
			nThreadsMin = atoi(tmp.c_str());

			// SECOND VALUE (If given)
			if (i + 1 < argc && argv[i + 1][0] != '-') {
				// Validate second number
				tmp = string(argv[++i]);
				if (!regex_match(tmp, uintNumRegex)) {
					cout << "Invalid nThread value format. Specify unsigned integer value or two if range." << endl;
					return 1;
				}

				// Parse C-str
				nThreadsMax = atoi(tmp.c_str());
			}
			else {
				nThreadsMax = nThreadsMin;
			}

		}
		else if (strcmp(argv[i], "-iso") == 0) {
			// ISO value is given
			if (i + 1 >= argc) {
				cout << "Missing iso value after -iso flag." << endl;
				return 1;
			}

			isoValueStr = string(argv[++i]);

			if (!regex_match(isoValueStr, anyNumRegex)) {
				cout << "Invalid iso value format. Please specify float." << endl;
				return 1;
			}
		}
		else if (strcmp(argv[i], "-iter") == 0) {
			// ISO value is given
			if (i + 1 >= argc) {
				cout << "Missing number of iterations after -iter flag." << endl;
				return 1;
			}

			string iterationsStr = string(argv[++i]);

			if (!regex_match(isoValueStr, uintNumRegex)) {
				cout << "Invalid iterations value format. Please specify unsigned integer." << endl;
				return 1;
			}

			iterations = atoi(iterationsStr.c_str());
		}
		else if (strcmp(argv[i], "-o") == 0) {
			// Output file path is specified
			if (i + 1 >= argc || argv[i + 1][0] == '-') {
				cout << "Missing file path after -o flag." << endl;
				return 1;
			}

			// Store the file name and offset iterator
			outFileStream = ofstream(argv[++i]);
			
			if (!outFileStream.good()) {
				cout << "Specified output file path is invaild." << endl;
			}
		}
		else if (strcmp(argv[i], "-data-type") == 0) {
			// Volume data type is specified
			if (i + 1 >= argc || argv[i + 1][0] == '-') {
				cout << "Missing type after -data-type flag." << endl;
				return 1;
			}

			// Data type is specified (char, uchar, short, ushort, int, uint, float, double)
			if (strcmp(argv[i + 1], "char") != 0 && strcmp(argv[i + 1], "uchar") != 0 && strcmp(argv[i + 1], "short") != 0 && strcmp(argv[i + 1], "ushort") != 0 && strcmp(argv[i + 1], "uint") != 0 && strcmp(argv[i + 1], "float") != 0 && strcmp(argv[i + 1], "double") != 0) {
				cout << "Invalid data type. Available data types: char, uchar, short, ushort, int, uint, float, double." << endl;
				return 1;
			}

			type = string(argv[++i]);
		}
		else {
			cout << "Unknown parameter: " << argv[i] << endl;
			return 1;
		}
	}

	if (inputFileStream.is_open() && (!customSizeSpecified || type.compare("") == 0 || isoValueStr.compare("") == 0)) {
		cout << "If custom volume is imported, you must input volume dimensions(-vol-dim), data type (-data-type) and iso value (-iso)." << endl;
		return 1;
	}
	// END PARAMETER PARSING
	

	if (benchmark) {
		if (type.compare("char") == 0) {
			benchmarkHandler(inputFileStream, outFileStream, size, voxSize, (char)atoi(isoValueStr.c_str()), nThreadsMin, nThreadsMax, iterations);
		}
		else if (type.compare("uchar") == 0) {
			benchmarkHandler(inputFileStream, outFileStream, size, voxSize, (unsigned char)atoi(isoValueStr.c_str()), nThreadsMin, nThreadsMax, iterations);
		}
		else if (type.compare("short") == 0) {
			benchmarkHandler(inputFileStream, outFileStream, size, voxSize, (short)atoi(isoValueStr.c_str()), nThreadsMin, nThreadsMax, iterations);
		}
		else if (type.compare("ushort") == 0) {
			benchmarkHandler(inputFileStream, outFileStream, size, voxSize, (unsigned short)atoi(isoValueStr.c_str()), nThreadsMin, nThreadsMax, iterations);
		}
		else if (type.compare("int") == 0) {
			benchmarkHandler(inputFileStream, outFileStream, size, voxSize, (int)atoi(isoValueStr.c_str()), nThreadsMin, nThreadsMax, iterations);
		}
		else if (type.compare("uint") == 0) {
			benchmarkHandler(inputFileStream, outFileStream, size, voxSize, (unsigned int)atoi(isoValueStr.c_str()), nThreadsMin, nThreadsMax, iterations);
		}
		else if (type.compare("float") == 0) {
			benchmarkHandler(inputFileStream, outFileStream, size, voxSize, (float)atoi(isoValueStr.c_str()), nThreadsMin, nThreadsMax, iterations);
		}
		else if (type.compare("double") == 0) {
			benchmarkHandler(inputFileStream, outFileStream, size, voxSize, (double)atoi(isoValueStr.c_str()), nThreadsMin, nThreadsMax, iterations);
		}
		else {
			benchmarkHandler(inputFileStream, outFileStream, size, voxSize, (float)atoi(isoValueStr.c_str()), nThreadsMin, nThreadsMax, iterations);
		}
	}
	else {
		if (!outFileStream.is_open()) {
			cout << "To extract the data the output file path is needed (-o)." << endl;
			return 1;
		}

		if (type.compare("char") == 0) {
			extractInputHandler(inputFileStream, outFileStream, size, voxSize, (char)atoi(isoValueStr.c_str()), nThreadsMax);
		}
		else if (type.compare("uchar") == 0) {
			extractInputHandler(inputFileStream, outFileStream, size, voxSize, (unsigned char)atoi(isoValueStr.c_str()), nThreadsMax);
		}
		else if (type.compare("short") == 0) {
			extractInputHandler(inputFileStream, outFileStream, size, voxSize, (short)atoi(isoValueStr.c_str()), nThreadsMax);
		}
		else if (type.compare("ushort") == 0) {
			extractInputHandler(inputFileStream, outFileStream, size, voxSize, (unsigned short)atoi(isoValueStr.c_str()), nThreadsMax);
		}
		else if (type.compare("int") == 0) {
			extractInputHandler(inputFileStream, outFileStream, size, voxSize, (int)atoi(isoValueStr.c_str()), nThreadsMax);
		}
		else if (type.compare("uint") == 0) {
			extractInputHandler(inputFileStream, outFileStream, size, voxSize, (unsigned int)atoi(isoValueStr.c_str()), nThreadsMax);
		}
		else if (type.compare("float") == 0) {
			extractInputHandler(inputFileStream, outFileStream, size, voxSize, (float)atoi(isoValueStr.c_str()), nThreadsMax);
		}
		else if (type.compare("double") == 0) {
			extractInputHandler(inputFileStream, outFileStream, size, voxSize, (double)atoi(isoValueStr.c_str()), nThreadsMax);
		}
		else {
			extractInputHandler(inputFileStream, outFileStream, size, voxSize, (float)atoi(isoValueStr.c_str()), nThreadsMax);
		}
	}
	
	return 0;
}

