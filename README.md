# Scattered Points Interpolation with Globally Smooth B-Spline Surface using Iterative Knot Insertion

The implementation of the paper "Scattered Points Interpolation with Globally Smooth B-Spline Surface using Iterative Knot Insertion". This code is tested on Windows and Linux, and it probably works on Mac.

## Build

Before build the code, you need to have CMake installed on your machine. To build the library and executable benchmark on Linux or macOS run:

```sh
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```
## Usage
You can run the examples:
### Benchmark models 
```bash
./Sparse_Interp_bin benchmarks delta per nbr outpath
```
where `delta` and `per` are two input parameters, ranging from 0 to 1, we recommand that `delta = 0.9` and `per = 0.5`. `nbr` is the number of sampled points, we suggest to choose a number no more than 500. `outpath` is the output path. For example, on Windows, you can run `./Sparse_Interp_bin benchmarks 0.9 0.5 100 C:\\`.


### Mesh models 
Similarly, you can get the 2 mesh model interpolation results in the paper by running

```bash
./Sparse_Interp_bin meshes modelname outpath
```
where `modelname` choose between `mask3kf.obj` and `tiger.obj`. For example, run  `./Sparse_Interp_bin meshes mask3kf.obj C:\\`.

### Reproduction of local non-smoothness
You can run the results of Figure 7 in the paper, by
```bash
./Sparse_Interp_bin reproduce outpath
```
It will reproduce the non-smoothness problem by inserting a lot of redundant knots into a fine surface, and reduce local energy by inserting knots to regions with large energy.

### Fixing local non-smoothness caused by too few data points
You can run the results of Figure 7 in the paper, by
```bash
./Sparse_Interp_bin local_energy outpath
```
It will produce the results in Figure 9.

## Outputs

After running the code, in `outpath` there will be a interpolation surface in `.obj` format, and for each interpolation surface, there is a corresponding `.csv` file recording its runtime, maximal interpolation error, and the number of control points. After running the benchmarks, the 6 benchmark models: Peak, Drop, Hyperbolic, Sinus, Bilinear and Snail will be output into `outpath`. After running the results in Figure 7 and Figure 9, the maximal energy for each iteration will also be written into the `.csv` files.

