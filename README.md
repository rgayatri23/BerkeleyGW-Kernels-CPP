# BerkeleyGW-kernels
The repo contains 2 compute kernels from the BerkeleyGW software package.
1. General Plasman Pole Self-Energy Simulation (GPP)
2. Full Frequency Self-Energy Simulation (FF)

The kernels in this repo are C++ ports or the orginal FORTRAN implementation in the BerkeleyGW software package.
The purpose of these kernels is to test the portability and performance of OpenMP and OpenACC programming models with C++ kernels.
Since there are multiple implementations of OpenMP and OpenACC provided by various compiler developers, the build system is tuned to accept parameters that define the compiler and implementation.

The kernels perform complex number arithemetic and use a CustomComplex class for complex number representation. This class is in the ComplexClass directory.
For Multi-Dimensional data structures, we define an arrayMD class that is available in the arrayMD directory.
These classes do not provide all the features of complex number arithemetic of MD arrays but contain only features that are used by the mini-apps.

## Build instructions
The repo consists of a CMake and a Makefile build system.

### CMake
The root directory consists of a CMakeLists.txt .
The CMake build can be used for Serial, OpenMP3.0 and OpenMP4.5 builds but not for OpenACC build.
The flags for OpenMP4.5 builds are only set for LLVM/Clang and IBM/XL compilers.
To build OpenMP4.5 version of the code with LLVM compiler use the follwing set of commands

```
  mkdir build && cd build
  cmake -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang -DOPENMP_TARGET=ON -DCMAKE_BUILD_TYPE=Release ../
```
By default the CMake build system will only build GPP, to build FF use `-DBUILD_FF=ON` .
For OpenMP3.0 version replace `OPENMP_TARGET` with `OPENMP` .

### Raw Makefile
There is a Makefile provided inside individual app directories which can be used by providing the following information
1. Compiler (comp=xl/clang/icc/pgi)
2. Programming model (OPENMP=y && OPENMP_TARGET=y or OPENACC=y)
3. Underlying architecture (GPU=y)
For example, to build OpenMP4.5 version of GPP with LLVM compiler use the follwing set of commands

```
  cd GPP/
  make COMP=clang OPENMP=y OPENMP_TARGET=y GPU=y
```
