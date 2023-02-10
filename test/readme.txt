Introduction
This repository includes executable programs and data for related algorithms in the paper "Hexagon-based Adaptive Hierarchies for Efficient Point-in-Spherical-Polygon Tests on GPUs". When the software and hardware environment meet requirements, you can run the bat file to verify the data in the article. The main bat files include testAHHG.bat and testGAHHG.bat, which are used to test the AHHG and the GAHHG method in the paper, respectively. Implementations of the other three comparison algorithms (GMAT, Lune(no prep), and Lune(prep)) are not included in this repository. They are provided by [1] and the programs are at the link https://github.com/ryanketzner/sphericalpolygon.


Directory
The source directory includes the source code files.
The test directory includes the test files.
The test/polygon directory includes polygon data for testing. The file format is customized. For details, see below.
The test/query directory contains the query point data for testing. The file format is csv.
The test/data_for_lune directory includes polygon files and query point files stored in the format used in [1] (note: only the data not provided by [1] is included).


Running requirments
The operating system of the program running environment is windows 10. And the operation of the GAHHG program requires a graphics card that supports CUDA 10 or higher.


Running method
Run testAHHG.bat to get the test result of the AHHG method. Performance statistics (including  preprocessing time, query time and space overhead, etc.) are recorded in the report.txt file. Similarly, run the testGAHHG.bat file to get the test results of the GAHHG method. But before starting the test, it checks to see if there is a GPU available. If not, the test is exited.


Polygon file format
The custom polygon file format is as follows.

POLYGON 1 [number of inner rings of the polygon]
0 [number of vertices of the outer ring]
[longitude] [latitude]
......
[longitude] [latitude]
[index of the inner ring] [number of vertices of the inner ring]
[longitude] [latitude]
......
[longitude] [latitude]
[index of the inner ring] [number of vertices of the inner ring]
[longitude] [latitude]
......

The value of longitude is [-180,180], and the value of latitude is [-90,90]. If a ring has n vertices, n+1 vertices are included in the file, where the first vertex is the same as the n+1 vertex. Here is an example file of a spherical triangle with an inner ring.

POLYGON 1 1
0 4
45.000023 30.000024
-45.000001 30.000001
0.000000 -35.000002
45.000023 30.000024
0 4
0.000000 -20.000001
-19.999999 19.999996
19.999999 20.000001
0.000000 -20.000001


Referrence
[1] Ketzner, R., Ravindra, V., and Bramble, M., 2022. A robust, fast, and accurate algorithm for
point in spherical polygon classification with applications in geoscience and remote sensing.
Computers and Geosciences, 167 (4), 615–624.