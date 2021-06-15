# Boundary-Sampled Halfspaces

![](figure/featured.jpeg)

[Xingyi Du](https://duxingyi-charles.github.io/), [Qingnan Zhou](https://research.adobe.com/person/qingnan-zhou/),  [Nathan Carr](https://research.adobe.com/person/nathan-carr/), [Tao Ju](https://www.cse.wustl.edu/~taoju/)
<br/>*ACM Transaction on Graphics (Proceedings of SIGGRAPH 2021)*<br/>

[`Project Page`](https://duxingyi-charles.github.io/publication/boundary-sampled-halfspaces/)

## Abstract

We present a novel representation of solid models for shape design. Like Constructive Solid Geometry (CSG), the solid shape is constructed from a set of halfspaces without the need for an explicit boundary structure. Instead of using Boolean expressions as in CSG, the shape is defined by sparsely placed samples on the boundary of each halfspace. This representation, called Boundary-Sampled Halfspaces (BSH), affords greater agility and expressiveness than CSG while simplifying the reverse engineering process. We discuss theoretical properties of the representation and present practical algorithms for boundary extraction and conversion from other representations. Our algorithms are demonstrated on both 2D and 3D examples.

## Code

- BSH_CLI, a command line program that extracts the boundary of the BSH shape defined by a set of halfspaces and sample points. 
- BSH_GUI, a GUI program for interactively solid modeling using our BSH representation.

The programs have been tested on macOS 11.4.

## Build

### Mac

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

The program `BSH_CLI` and `BSH_GUI` will be generated in the `build` subdirectory.

## Usage

### BSH_CLI

Usage: 

    ./BSH_CLI [OPTIONS] config_file output_grid_file

Positionals:
- config_file (REQUIRED): Configuration file, specifying input halfspaces and samples. See `/example/xxx/input/config.json` for examples.
- output_grid_file (REQUIRED): Output grid file, containing the extracted BSH boundary and other related information. See `/example/xxx/output/result.grid` for examples.

Options:
- -h,--help : Print help message and exit.
- -G,--grid-file (REQUIRED): Grid spec file, specifying bounding box and grid resolution for Marching Cube.
- -P,--param-file (REQUIRED): Parameter spec file, specifying parameters for BSH boundary extraction algorithm. For default parameters, see `param.json` files under `/example`.
- -A,--arr-algo (REQUIRED): Arrangement algorithm. Currently, only `mesh` ([this method](https://github.com/gcherchi/FastAndRobustMeshArrangements)) is supported.

Example:

    ./BSH_CLI ../examples/figure16/tori/input/grid_64.json -P ../examples/figure16/tori/input/param.json -A mesh ../examples/figure16/tori/input/config.json  ../examples/figure16/tori/output/result.grid


## Examples

`/example` directory contains input data and output results used to generate the corresponding figures in the BSH paper.