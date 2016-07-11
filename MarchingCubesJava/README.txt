USAGE:
This script may be executed in either benchmark or extract mode. Mode is specified by the first parameter [benchmark, extract].
Parameters:
    -input-vol       Specifies path to the input volume. If this parameter is set volume dimensions(-vol-dim), data type(-data-type) and iso value(-iso) must also be given.
    -vol-dim         Specifies the generated/read volume dimensions. Dimensions should be given as unsigned integers in format; -vol-dim X Y Z.
    -data-type       Specifies the input file or generated data type. Options [char, uchar, short, ushort, int, uint, float, double].
    -vox-dim         Specifies voxel dimensions used in mesh construction. Dimensions should be given as floating point numbers in format: -vox-dim X Y Z.
    -nThread         Number of threads used in Marching cubes algorithm.This parameter can be either given as a single unsigned integer value or two unsigned integer values in benchmark mode, specifying the range of thread executions that will be tested.
    -iter            Used only in benchmark mode to determine how many iterations should be executed for each configuration.
    -iso             Isovalue that is used as a threshold for determining active voxels. Type should match the data type.
    -o               Path to output file. In extract mode the mesh is written to file in .obj format [required]. In benchmark mode the results are written to file.
