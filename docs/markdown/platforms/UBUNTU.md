## Ubuntu
It is possible to build GridPACK on Ubuntu by following instructions for
[CentOS/RHEL 7](CENTOS.md) or the instructions for the [basic
build](../BASIC_BUILD.md). However, it is also possible to reduce many of the
steps in the installation process using the package manager.
The Global Arrays library must still be built by hand but MPI, Boost and PETSc
can all be installed using the package manager. This
example was performed on a [Virtual Box](https://www.virtualbox.org/) instance
running a clean install of [Ubuntu Linux
18.04](http://releases.ubuntu.com/18.04/) (LTS). Virtual Box is not required
and this build should work on Linux systems
using the Ubuntu operating system. **If you are installing packages, you will need
super user or sudo privileges for this installation**. You are most likely to
need to do the package installation if you are working in a newly created VM or
are configuring a workstation from scratch. If you are using a workstation that
is already configured for scientific computations, many of these packages may
already be installed.

A few basic development packages will be necessary. Install them with

```
sudo apt-get install git build-essential devscripts equivs
```


## Install Required Packages
Before installing the remaining packages, you **must** first download GridPACK
onto your computer and put it in a convenient directory. Information on how to
download and install GridPACK can be found
[here](../required/GRIDPACK.md#obtaining-gridpack-from-source).

All the [prerequisite packages](../required/PREREQUISITES.md)
can be installed using the packaging
information in the source. Run this command in the top directory of the cloned
repository:

```
sudo mk-build-deps -i
```

This will make and install a virtual package called
`gridpack-build-deps` with most of the dependencies needed to build
GridPACK.

## Build and Install Global Arrays

Get the [Global Arrays](http://hpc.pnl.gov/globalarrays/) source from the
[release page](https://github.com/GlobalArrays/ga/releases/tag/v5.8.2). We
recommend using a version greater than 5.6. Unpack GA
in a convenient location. In the unpacked directory, configure, build, and
install using

```
    ./configure --enable-i4 --without-blas --enable-cxx --with-mpi-ts --disable-f77 \
        --prefix="/home/gridpack/software/ga-5.8"
    make
    make install 
```

Be sure and replace the install directory
with an appropriate directory on your local system.

## Build and Test GridPACK 

[Obtain the GridPACK release or development
code](../required/GRIDPACK#obtaining-gridpack-source) and put
it in a convenient directory.

From the top-most GridPACK source directory, do the following:

```
  mkdir src/build
  cd src/build
```

Modify this configuration recipe to fit your needs. The main changes you need to
make are to modify the directories to reflect your local file structure.

```
  rm -rf CMake*
  export CC=gcc
  export CXX=g++
  
  cmake  \
    -D PETSC_DIR:STRING="/usr/lib/petsc" \
    -D PARMETIS_DIR:PATH="/usr" \
    -D GA_DIR:STRING="/home/gridpack/software/ga-5.8" \
    -D MPI_CXX_COMPILER:STRING="mpicxx" \
    -D MPI_C_COMPILER:STRING="mpicc" \
    -D MPIEXEC:STRING="mpiexec" \
    -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
    -D GRIDPACK_TEST_TIMEOUT:STRING="20" \
    -D USE_GLPK:BOOL=ON \
    -D GLPK_ROOT_DIR:PATH="/usr" \
    -D BUILD_SHARED_LIBS:BOOL=OFF \
    -D CMAKE_INSTALL_PREFIX:PATH="/home/gridpack/gridpack-3.3/src/build/install" \
     CC=gcc CXX=g++ ..
  
  make
  make install
```

Again, be sure and modify the Global Array directory and the GridPACK
installation directories to match the
appropriate locations on your local system.

If compilation is successful, the [unit
tests and/or example applications](../required/GRIDPACK.md#running-tests)
can be run.

## Removing Required Software

If no longer needed, the required packages can be removed using

```
sudo apt-get purge gridpack-build-deps
sudo apt autoremove
```

