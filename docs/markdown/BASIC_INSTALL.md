# GridPACK installation instructions
This document provides a step-by-step guide to install GridPACK and its
dependencies. They are meant to be used with GridPACK's develop branch. The
installation instructions have been tested on Linux, MacOS, and Ubuntu. Before
building GridPACK using these instructions, you will need to make sure that
CMake is available on your system (version newer than 3.5.0) and that you have
a GNU compiler installed. You will also need to have MPI installed on your
platform. Most clusters already have MPI available, but if you are building
GridPACK on a workstation or a virtual machine, you may need to install or build
MPI on your own. More information on configuring your Linux platform
can be found [here](required/LINUX_BASICS.md).

The example below can be adapted to build GridPACK on many platforms. More
specific examples, with detailed instructions, for builds using different
variants of Linux can be found on the links below.

* [CentOS7 or RHEL 7](platforms/CENTOS7.md)
* [Ubuntu](platforms/UBUNTU.md)
* [Ubuntu PPA Install](platforms/UBUNTU_PPA.md)
* [Debian](platforms/Debian.md)
* [Linux Cluster](platforms/RC_CLUSTER.md)
* [Notes](../notes): GridPACK has been built on several leadership platforms.
  These builds have typically been fairly complicated and require assistance
  from the support staff. Notes on these builds for individual platforms can be
  found in this directory and may be of use to those attempting to build on
  these platforms or similar platforms.

## Submodule(s)

If building a verions of GridPACK that has been cloned directly from the Github
repository, it is first necessary to download some submodules. This can be done
using the command

```
git submodule update --init
```

The CMake configuration will fail if this command has not been run. If you are
building GridPACK from one of the releae tarballs, the submodules will already
be included in the tarball and this command is not necessary.

## Boost 1.78.0
Step 1. Download boost 1.78.0
```
wget https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.gz
```
Step 2. Untar file
```
tar -xvf boost_1_78_0.tar.gz
```
Step 3. `cd boost_1_78_0`

Step 4. Build Boost with the needed libraries and set the installation directory. Here, we'll use `install_for_gridpack` as the installation directory.

```
./bootstrap.sh --prefix=install_for_gridpack --with-libraries=mpi,serialization,random,filesystem,system

```
Step 5. Edit `project-config.jam` file and add a line to it with the text
`using mpi ;`

Step 6. Install Boost. This is done by the following two commands
```
./b2 -a -d+2 link=static stage
```
and
```
./b2 -a -d+2 link=static install
```

**Note:** To build shared libraries for Boost, use the flag `link=shared` instead of `link=static`.

Step 7. If boost installation goes through correctly, then you should see `include` and `lib` subdirectories in your installation directory `install_for_gridpack`

If you run into difficulties more information on building and installing Boost
can be found [here](required/BOOST.md)

## Install Global arrays
Step 1. Download GA-5.8
```
wget https://github.com/GlobalArrays/ga/releases/download/v5.8/ga-5.8.tar.gz
```
Step 2. Untar file
```
tar -xvf ga-5.8.tar.gz
```
Step 3. `cd ga-5.8`

Step 4. Configure ga and set the installation directory
```
./configure --with-mpi-ts --disable-f77 --without-blas --enable-cxx --enable-i4 --prefix=${PWD}/install_for_gridpack
```
**Note:** To build shared libraries for GA, add the flag `--enable-shared`

Step 5. Compile GA and install
```
make -j 10 install
```
Step 6. If the compilation succeeds, `include`, `lib`, and `bin` directories should be created in your installation directory `${PWD}/install_for_gridpack`

More information on building and installing GA can be found
[here](required/GLOBAL_ARRAYS.md)

## Install PETSc 3.16.4
Step 1. Download PETSc release

```
git clone https://gitlab.com/petsc/petsc.githttps://gitlab.com/petsc/petsc.git
git checkout v3.16.4
```
Step 2. `cd petsc`

Step 3. Set environment variables needed for PETSc
```
export PETSC_DIR=${PWD}
export PETSC_ARCH=build-dir
```
Step 4. Build PETSc. GridPACK needs a few other libraries installed with PETSc. We'll install them during configure stage of PETSc.
```
./configure --download-superlu_dist --download-metis --download-parmetis --download-suitesparse --download-f2cblaslapack --download-cmake --prefix=${PWD}/install_for_gridpack
```
**Note:** To build shared libraries for PETSc, add the flag `--with-shared-libraries=1`

**Note:** `--download-f2cblaslapack` and `--download-cmake` is not required if the system has a functional BLAS/LAPACK compatible with PETSc, and the `CMake` version is 3.18.1 or higher. On EIOC-Ubuntu8 you'll need to add both these configuration options.

This step will take a few minutes.

Step 5. Once the configuration is complete, run the following commands to compile and test PETSc
```
make
make install
make check
```

Step 6. Once PETSc installation is complete and all tests pass after running `make check`, you'll see three directories `include`, `lib`, and `bin` under the installation directory `${PWD}/install_for_gridpack`

More information on building and installing PETSc can be found
[here](required/PETSC.md). In the unlikely event that PETSc cannot download and
build ParMETIS, more information of building this library can be found
[here](required/PARMETIS.md).

## Building GridPACK
Now that we have installed the needed dependencies, we can proceed with downloading and installing GridPACK.

Step 1. Download GridPACK
```
git clone https://github.com/GridOPTICS/GridPACK.git
```
Step 2. Change branch to develop
```
git checkout develop
```
Step 3. `cd GridPACK/src`

Step 4. Create build directory `mkdir build`, and cd to it - `cd build`

Step 5. Copy the following GridPACK configuration commands to a file (name it for e.g. `build.sh`)
```
rm -rf CMake*                                                                               
cmake \                                                                                     
   -D GA_DIR:STRING="/home/abhy245/software/ga-5.8/install_for_gridpack" \                
   -D BOOST_ROOT:STRING="/home/abhy245/software/boost_1_78_0/install_for_gridpack" \    
     -D Boost_DIR:STRING="/home/abhy245/software/boost_1_78_0/install_for_gridpack/lib/cmake/Boost-1.78.0" \   
   -D BOOST_LIBRARYDIR:STRING="/home/abhy245/software/boost_1_78_0/build/lib"\   
   -D PETSC_DIR:PATH="/home/abhy245/software/petsc/install_for_gridpack" \                  
   -D MPI_CXX_COMPILER:STRING='mpicxx' \     
   -D MPI_C_COMPILER:STRING='mpicc' \                                         
   -D MPIEXEC:STRING='mpiexec' \                                                            
   -D GRIDPACK_TEST_TIMEOUT:STRING=30 \                                                     
   -D CMAKE_INSTALL_PREFIX:PATH="/home/abhy245/software/GridPACK/src/install" \             
   -D CMAKE_BUILD_TYPE:STRING=Debug \                                                       
    ..   
```
**Note:** To build shared libraries, add `-D BUILD_SHARED_LIBS=YES` to the abovve commands.

**Note:** Change the paths for `GA_DIR`,`BOOST_DIR`,`PETSC_DIR` to point to the installation directories for `GA`, `Boost`, and `PETSc`, respectively.

**Note:** Change the path for `CMAKE_INSTALL_PREFIX` to your choice. This is the directory where GridPACK will be installed.

**Note:** If you do not use `--prefix=${PWD}/install_for_gridpack` flag for PETSc install then PETSc will write the installation files in `$PETSC_DIR/$PETSC_ARCH`. So, you'll need to add `-D PETSC_ARCH=<PETSC_ARCH_NAME>` to the above CMake commands.

Step 6. Make the file executable `chmod +x build.sh`

Step 7. Run `build.sh`
```
./build.sh
```
Step 8. Compile and install
```
make -j 10 install
```
During installation, you may see several warnings generated. This is typical of GridPACK install. You can ignore them.

Step 9. Test GridPACK install by running tests.
```
make test
```
You may perhaps see several of the tests failing. This is fine. We'll test in the next step whether the dynamics simulation application is running correctly.

Step 10. Test GridPACK install by running dynamics simulation
```
cd applications/dynamic_simulation_full_y
cp input_145.xml input.xml
./dsf.x
```

More information on configuring and building GridPACK can be found
[here](required/GRIDPACK.md)


