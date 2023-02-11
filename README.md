# GAHHG

## Introduction

This repository contains the souce code, executable files and data for the article "Hexagon-based Adaptive Hierarchies for Efficient Point-in-Spherical-Polygon Tests on GPUs". When the software and hardware conditions are met, readers can verify the data in the article by compiling source code or directly running bat files. The main bat file is testGAHHG.bat, which is used to test the GAHHG algorithm in this article. The performance of the comparison algorithm AHHG in the paper can be tested by running testAHHG.bat. For the other three comparison algorithms (GMAT, Lune(no prep), and Lune(prep)), they are provided by [1] and the relating data and programs are at the link https://github.com/ryanketzner/sphericalpolygon. For easy comparison, this repository includes some polygon files and query point files stored in the format they use.


## Directory
```
├── source
│   ├── scan
│   └── work 
└── test
    ├── polygon
    ├── query
    └── data_for_lune
```

The `source` directory contains the source code for the GAHHG implementation, which is a Visual Studio 2019 project that uses the CUDA 11.6 runtime. 

The `source/scan` directory contains the parallel scan source code used by the GAHHG algorithm, which comes from https://github.com/mattdean1/cuda.

The `source/work` is the working directory of the VS project.

The `test` directory includes compiled executable files for AHHG and GAHHG algorithms, as well as batch files and data for testing.

The `test/polygon` directory contains polygon files for testing. The polygon file format is customized. For details, see below.

The `test/query` directory contains query point files for testing. The query point file format is csv.

The `test/data_for_lune` directory includes polygon files and query point files stored in the format used in [1] (note: only the data not provided by [1] is included).

## Build the project

To build the project, it requires the Windows 10 operating system, Visual Studio 2019 and CUDA Toolkit 11.6.

## Reproduce the result in the paper

To reproduce the result, it requires the Windows 10 operating system and a graphics card that supports CUDA 11 or higher. Run testGAHHG.bat file in the test directory to get the test result of the GAHHG method. Performance statistics (including  preprocessing time, query time and space overhead, etc.) are recorded in the report.txt file. Similarly, run the testAHHG.bat file in the test directory to get the test results of the AHHG method. Before starting the GAHHG test, it checks to see if there is a GPU available. If not, the test exits.

## Polygon file format

The custom polygon file format is as follows.

POLYGON 1 [number of inner rings of the polygon]

0 [number of vertices of the outer ring + 1]

[longitude] [latitude]

......

[longitude] [latitude]

[index of the inner ring] [number of vertices of the inner ring + 1]

[longitude] [latitude]

......

[longitude] [latitude]

[index of the inner ring] [number of vertices of the inner ring + 1]

[longitude] [latitude]

......

The longitude range is -180° ~ 180° and the latitude range is -90° ~ 90°. If a ring has n vertices, n+1 vertices are included in the file, where the first vertex is the same as the (n+1)th vertex. Below is an example file of a spherical triangle with an inner ring.

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

## Referrence

[1] Ketzner, R., Ravindra, V., and Bramble, M., 2022. A robust, fast, and accurate algorithm for
point in spherical polygon classification with applications in geoscience and remote sensing.
Computers and Geosciences, 167 (4), 615–624.
