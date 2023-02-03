Building GridPACK is relatively straightforward on [Debian 9
(stretch)](https://www.debian.org/) systems.  At the time of writing,
[Debian 9](https://www.debian.org/) was the current stable distribution.
There is no need to build any
[prerequisite software](../required/PREREQUISITES.md).  All software can be
installed from Debian package repositories.  This documents a GridPACK build on
Debian 9 (stretch) installed on a
[VirtualBox](https://www.virtualbox.org/) virtual machine (VM) using the
[complete installation image](https://www.debian.org/CD/). Virtual Box is not
required and this build should work on Linux systems using the Debian operating
system.

## System Preparation

**You will need super user or sudo privileges for this installation**. You will
not be able to edit files
`/etc/apt/sources.list` or use utilities such as `apt-get`
without them.  Here, the `sudo` command is used to perform super user
activities. This command was installed automatically in this case. Your system
may be different as it was in
[at least one other case](https://github.com/GridOPTICS/GridPACK/issues/27).

Starting with a **clean** installation, add `contrib` and
`non-free` components of the Debian distribution to apt sources (ParMETIS
is in non-free). Edit `/etc/apt/sources.list` and make the main
repository line look like this:

```
 deb http://ftp.us.debian.org/debian/ stretch main contrib non-free
 deb-src http://ftp.us.debian.org/debian/ stretch main contrib non-free
```

It may be necessary to change the file permissions before editing:

```
 sudo chmod +w /etc/apt/sources.list
```

Refresh the system package lists with

```
 sudo apt-get update
```

## Prerequisite Installation

###General

Install a C++ compiler, [CMake](https://cmake.org/), and
[Git](https://git-scm.com/):

```
 sudo apt-get install git cmake g++
```

### Boost

Install necessary [Boost](http://www.boost.org/) libraries:

```
 sudo apt-get install libboost-dev libboost-mpi-dev \
    libboost-random-dev libboost-filesystem-dev libboost-system-dev
```

This will also install the default MPI implementation (OpenMPI), including
compiler wrappers.

### PETSc

Install the real-valued version of [PETSc](https://www.mcs.anl.gov/petsc/) with

```
 sudo apt-get install petsc-dev
```

or the complex-valued version with

```
 sudo apt-get install libpetsc3.7.5-dev
```

### Global Arrays

[Global Arrays](https://github.com/GlobalArrays/ga) has lots of dependencies.
Unfortunately, installing the Debian GA package (version
`5.4~beta~r10636+dfsg-5`) does not enforce any of them. Installing
[https://www.mcs.anl.gov/petsc/ PETSc] first will install most of them. This
should complete the GA installation:

```
 sudo apt-get install libglobalarrays-dev libarmci-mpi-dev
```

### ParMETIS

Install [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
using

```
 sudo apt-get install libparmetis-dev libmetis-dev
```

### GNU Linear Programming Kit

[GLPK](https://www.gnu.org/software/glpk/) is optional and can be installed with

```
 sudo apt-get install libglpk-dev
```

## GridPACK Configuration and Build

[Obtain the GridPACK release or development code](https://github.com/GridOPTICS/GridPACK/releases)
and put it in a convenient directory, like `$HOME/gridpack/src`. The top level
GridPACK directory is denoted below by the variable `$HOME`.

It is a good idea to build GridPACK in a separate directory under the GridPACK
source tree. The example below assumes that a directory called `build`
has been created under `$GRIDPACK/src` and that you have cd'd into this
directory:

```
  cd $GRIDPACK/src
  mkdir build
  cd build
```

Configure GridPACK as follows

```
CC=gcc
CXX=g++
CFLAGS=-pthread
CXXFLAGS=-pthread
export CC CXX CFLAGS CXXFLAGS

cmake \
    -D PETSC_DIR:STRING="/usr/lib/petsc" \
    -D PARMETIS_DIR:PATH="/usr" \
    -D GA_EXTRA_LIBS:STRING="-lscalapack-openmpi -lblacs-openmpi -llapack -lblas -lgfortran" \
    -D MPI_CXX_COMPILER:STRING="mpicxx" \
    -D MPI_C_COMPILER:STRING="mpicc" \
    -D MPIEXEC:STRING="mpiexec" \
    -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
    -D GRIDPACK_TEST_TIMEOUT:STRING=30 \
    -D USE_GLPK:BOOL=ON \
    -D GLPK_ROOT_DIR:PATH="/usr" \
    -D CMAKE_BUILD_TYPE:STRING=Release \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
    -D CFLAGS="-pthread" FCFLAGS="-pthread" CXXFLAGS="-pthread" \
    ..
```

Then, build 

```
   make
```

If compilation is successful, the
[GriPACK test suite and/or examples](../required/GRIDPACK.md#running-tests)
can be run.
