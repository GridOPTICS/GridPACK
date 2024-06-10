# Quick guide for GridPACK-EMT
## Installation
GridPACK depends on three external packages (Boost, Global Arrays (GA), and PETSc). This section guides through the installation of these packages via a script (`install_gridpack_deps.sh`) we've developed.

### Prerequisites
Before beginning with the GridPACK installation process, we need to have some prerequisite libraries installed. These libraries are the fundamental building blocks of parallel applications and their build process.  

- Make sure you have a functional MPI installation. You can use either OpenMPI or MPICH. MPI is a prerequisite for GridPACK and its dependencies.
- Install CMake ver. 3.18.1 or higher. This is required for PETSc dependencies. You could also have PETSc install CMake. This can be done by adding `--download-cmake` to the PETSc configure. (line 118 in install_gridpack_deps.sh).
- You also need to have functional `Blas/Lapack` libraries that can be found by the installer. If you want PETSc to install Blas/lapack then add `--download-f2cblaslapack=1` to the PETSc configure line (line 118 in install_gridpack_deps.sh).

### Install GridPACK dependencies
1. Go to the root GridPACK directory.
```
cd $GRIDPACK_ROOT_DIR
```
Here, `$GRIDPACK_ROOT_DIR` is the root GridPACK directory


3. Assuming you have installed MPI and CMake, source `install_gridpack_deps.sh` to install all GridPACK dependencies. This will install all the GridPACK dependencies Boost (ver. 1.81), Global Arrays (ver. 5.8), and PETSc (ver. 3.20.1) in `$GRIDPACK_ROOT_DIR/external-dependencies`.
```
source install_gridpack_deps.sh
```

If you wish to install the libraries at another location then set `GP_EXT_DEPS` environment variable to the location you want to install the libraries. For example,
```
export GP_EXT_DEPS=/usr/local
```
will install all the external packages in `usr/local`

4. Once the dependencies get installed, it is best to set some environmental variables and update library paths so that the linker is always able to find the libraries (this is only to do if you have installed the dependencies in a non-standard location)
```
# Root GridPACK directory                                                                                        
#export GRIDPACK_ROOT_DIR=<location_of_GridPACK_directory>                      
# GridPACK external dependencies                                                                            
#export GP_EXT_DEPS=${GRIDPACK_ROOT_DIR}/external-dependencies                                                     
# GridPACK install directory                                                                                       
#export GRIDPACK_DIR=${GRIDPACK_ROOT_DIR}/src/install                                                              
# update LD_LIBRARY_PATH so that boost,ga, and petsc are on it                                                     
#export LD_LIBRARY_PATH=${GP_EXT_DEPS}/boost_1_78_0/install_for_gridpack/lib:${GP_EXT_DEPS}/ga-5.8/install_for_gri\
dpack/lib:${GP_EXT_DEPS}/petsc/install_for_gridpack/lib:${LD_LIBRARY_PATH}
```
### Installation errors
Here are some of the installation errors we've seen on different systems.
1. <b>PETSc cannot find Blas/Lapack:</b> This is because PETSc cannot find a functional Blas/Lapack. The error looks like this...
```
--------------------------------------------------------------------------------------------
  Could not find a functional BLAS. Run with --with-blas-lib=<lib> to indicate the library
  containing BLAS.
  Or --download-f2cblaslapack=1 to have one automatically downloaded and installed
*********************************************************************************************

Building PETSc 3.20.1
makefile:24: /home/abhy245/software/GridPACK/external-dependencies/petsc/build-dir/lib/petsc/conf/petscrules: No such file or directory
make[1]: *** No rule to make target '/home/abhy245/software/GridPACK/external-dependencies/petsc/build-dir/lib/petsc/conf/petscrules'.  Stop.
...
...
```
To fix this, add `--download-f2cblaslapack=1` to the PETSc configure line (line 118 in install_gridpack_deps.sh) and re-run the installation script.

2. <b>PETSc cannot find CUDA libraries</b>

If your system has NVIDIA CUDA libraries then PETSc may throw the following warning when completing its installation. This is a harmless warning. <b> You can ignore it. </b>

```
...
Running PETSc check examples to verify correct installation
Using PETSC_DIR=/people/abhy245/software/GridPACK/external-dependencies/petsc and PETSC_ARCH=build-dir
Possible error running C/C++ src/snes/tutorials/ex19 with 1 MPI process
See https://petsc.org/release/faq/
--------------------------------------------------------------------------
The library attempted to open the following supporting CUDA libraries,
but each of them failed.  CUDA-aware support is disabled.
libcuda.so.1: cannot open shared object file: No such file or directory
libcuda.dylib: cannot open shared object file: No such file or directory
/usr/lib64/libcuda.so.1: cannot open shared object file: No such file or directory
/usr/lib64/libcuda.dylib: cannot open shared object file: No such file or directory
If you are not interested in CUDA-aware support, then run with
--mca mpi_cuda_support 0 to suppress this message.  If you are interested
in CUDA-aware support, then try setting LD_LIBRARY_PATH to the location
of libcuda.so.1 to get passed this issue.
--------------------------------------------------------------------------
lid velocity = 0.0016, prandtl # = 1., grashof # = 1.
Number of SNES iterations = 2
Possible error running C/C++ src/snes/tutorials/ex19 with 2 MPI processes
See https://petsc.org/release/faq/
--------------------------------------------------------------------------
```

<b>Note:</b> If the error is only with PETSc install then you can have the installation script just install PETSc. This can be done by setting the `install_boost` flag (line 7) and `install_ga` flag (line 8) to `false`. This way, only PETSc installation will be executed.

## Installing GridPACK
Once the dependencies are installed, you can install the GridPACK library by running
```
source install_gridpack.sh
```
This will build and install the GridPACK library. The build files will be in `$GRIDPACK_ROOT_DIR/src/build` and the install files will be in `$GRIDPACK_ROOT_DIR/src/install`

## Executing GridPACK-EMT
To be completed.

