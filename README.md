# sparse data interpolation

The implementation of the paper "Scattered Points Interpolation with Globally Smooth B-Spline Surface using Iterative Knot Insertion".

## Build

To build the library and executable benchmark on Linux or macOS run:
```sh
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```
## Usage
You can run the examples:
- To run the benchmark models, you can 
```bash
./Sparse_Interp_bin benchmarks delta per nbr outpath
```
