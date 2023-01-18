## Linux Cluster with Infiniband

This page describes a GridPACK build on a Linux cluster using an Infiniband network
for communication. This build is similar to the one described for building
GridPACK on a Red Hat Enterprise Linux (RHEL) workstation, so many of the
options will be similar. The scripts described below should work for a variety
of RHEL systems, but we can verify that they worked on a Linux cluster running
RHEL 6.6 with OpenMPI 1.8.3 and using the GNU 4.9.2 compilers. The specific
versions of the libraries used to build GridPACK are Boost 1.65.0, Global Arrays
5.6.3 and PETSc 3.7.6. Other versions will likely work with little or no
modification, but these versions have actually been built using the scripts
below. We have also used these scripts to build GridPACK using Intel compilers
and intelmpi. A separate Boost configuration script is needed, but apart from
modifying directory locations and the name of the `PETSC_ARCH` variable,
the same configuration scripts can be used for PETSc, Global Arrays and
GridPACK. We have verified that the Intel build works on the same versions of
all libraries used in the GNU build with the Intel 17.0.4 compilers and
2018.0.128 version of intelmpi.

This build is very similar to the build for a [[Building_on_RHEL | workstation
using the RHEL operating system]]. The main difference is that this build puts
all libraries in separate directories after building, while the workstation
build puts all files in a single common directory. If you prefer to do you
builds so that all libraries and include files end up in a single directory, you
should review the workstation build for details on how to do this.

## Modules

If your system supports modules, then you can use these to configure your
environment so that the software needed to build GridPACK is in your path. For
building GridPACK using
GNU compilers and the OpenMPI version of MPI, use the modules

```
    module purge
    module load gcc/4.9.2
    module load openmpi/1.8.3
    module load python/2.7.3
    module load cmake/3.8.2
```

The details of these modules may vary on your system. In addition, set the
environment variables

```
    setenv CC gcc
    setenv CFLAGS "-pthread"
    setenv CXX g++
    setenv CXXFLAGS "-pthread"
    setenv FC gfortran
    setenv FFLAGS "-pthread"
```

The PETSc build will complain about these environment variables being set, but
it will just override them and use its own settings.

If you are interested in an Intel build, then try loading

```
    module purge
    module intel
    module intelmpi
    module python/2.7.3
    module cmake
```

and set the environment variables

```
    setenv CC icc
    setenv CFLAGS "-pthread"
    setenv CXX icpc
    setenv CXXFLAGS "-pthread"
    setenv FC ifort
    setenv FFLAGS "-pthread"
```

After setting modules, the output of the `module list` command for
the GNU environment should be

```
    Currently Loaded Modulefiles:
      1) gcc/4.9.2                3) python/2.7.3
      2) openmpi/1.8.3(default)   4) cmake/3.8.2
```

You will get something similar for Intel. Modules should guarantee that the
compilers and MPI libraries are all consistent. If you are not using modules,
then you will need to work with your system administrator to make sure that your
environment has a consistent set of compilers and MPI libraries.

## Boost

A Boost tar file can be obtained from the [download
page](https://www.boost.org/users/download/).
Information on downloading files can be obtained
[here](REQUIRED_SOFTWARE.md#linux-basics). Boost can configured and built
using the GNU environment as follows:

```
    echo "using mpi ;" > ~/user-config.jam
    sh ./bootstrap.sh \
        --prefix="/pic/projects/gridpack/software_new/boost_1_65_0</span>" \
        --without-icu \
        --with-toolset=gcc \
        --without-libraries=python
    ./b2 -a -d+2 link=static stage
    ./b2 -a -d+2 link=static install
    rm ~/user-config.jam
```

The argument to `--prefix` will need to be modified to reflect the
true location of the Boost libraries on your system. Make sure you include all
spaces in the above script exactly as written.

For Intel compilers and intelmpi, the script needs to be modified to

```
    echo "using mpi : mpicxx ;" > ~/user-config.jam
    sh ./bootstrap.sh \
        --prefix="/pic/projects/gridpack/software_intel/boost_1_65_0</span>" \
        --without-icu \
        --with-toolset=intel-linux \
        --without-libraries=python
    ./b2 -a -d+2 link=static stage
    ./b2 -a -d+2 link=static install
    rm ~/user-config.jam
```

Note that the argument to `--with-toolset` is changed and the
`mpicxx` wrapper has been added to the first line (this needs to be
done because intelmpi does not supply a compiler wrapper call `mpic++`).
Again, make sure you include all spaces in the above script exactly as written.

## PETSc

Tar files containing the complete PETSc library can be obtained from the
[PETSc download page](https://www.mcs.anl.gov/petsc/download/index.html).

In this example, PETSc was configured to include
[SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/). This requires that
ParMETIS and METIS be included also, so these libraries do not have to built
separately. The configuration script should be run in the top-level PETSc
directory. Once you have configure PETSc, there will be a directory
`linux-openmpi-gnu-cxx-complex-opt` that will contain all the PETSc
include and library files. Note that the value of the `PETSC_ARCH`
variable used to configure GridPACK (below) ''must'' match the value of
`PETSC_ARCH` used to configure PETSc.

```
    python ./config/configure.py \
      PETSC_ARCH=linux-openmpi-gnu-cxx-complex-opt</span> \
      --with-scalar-type=complex \
      --with-fortran-kernels=generic \
      --download-superlu_dist \
      --download-superlu \
      --download-parmetis \
      --download-metis \
      --download-suitesparse \
      --download-fblaslapack \
      --with-clanguage=c++ \
      --with-shared-libraries=0 \
      --with-x=0 \
      --with-mpiexec=mpiexec \
      --with-debugging=0
```

When you run this script, the PETSc build process will prompt you on what to do
next. Just copy and paste these instructions into the Linux command prompt. The
first instruction after configuring is

```
    xxx=========================================================================xxx
     Configure stage complete. Now build PETSc libraries with (gnumake build):
       make PETSC_DIR=/pic/projects/gridpack/software_new/petsc-3.7.6
PETSC_ARCH=linux-openmpi-gnu-cxx-complex-opt all
    xxx=========================================================================xxx
```

Follow the instructions and copy

```
    make PETSC_DIR=/pic/projects/gridpack/software_new/petsc-3.7.6
PETSC_ARCH=linux-openmpi-gnu-cxx-complex-opt all
```

into the command prompt. Additional instructions will appear at the end of the
make. If the make command succeeds, you can stop but you may want to run some of
the tests to insure that the build is correct.

If you are using Intel compilers and intelmpi, you should modify the name of
`PETSC_ARCH` variable to reflect this (although the build will still work
if you don't).

If you run into problems configuring or building PETSc, some additional tips can
be found [[Troubleshooting GridPACK Builds#PETSc | here]].

## Global Arrays

Download GA into a directory from the
[GA download page](https://github.com/GlobalArrays/ga/releases).
Since it is possible to configure GA in different
ways, it is often a good idea to create a build directory below the top level GA
directory to configure a build of GA. In this example, this directory is
`build_pr`. After creating the directory, cd into it and configure with

```
    cd $(GA_HOME)/build_pr
    ../configure --enable-i4 --enable-cxx --with-mpi-pr --without-blas
--disable-f77 --prefix=/pic/projects/gridpack/software_new/ga-5-6/build_pr
    make
    make install
```

This script assumes that there is a build directory called build_pr that is
beneath GA home directory. Be sure and include the `--without-blas`
directive, otherwise the configure will pick up the BLAS library from the
environment and this will conflict with BLAS that is downloaded and built by
PETSc. This script is configuring GA to use the progress ranks runtime, which is
recommended for large scale simulations. The `USE_PROGRESS_RANKS`
variable must be set to `TRUE` when configuring GridPACK (see below).
Also, because the progress ranks runtime requisitions 1 MPI process per SMP node
to manage computations, the actual number of processes available to the
calculation is less than the number requested when calling `mpirun`.
Because of this you must use at least 2 processes to run GridPACK applications
with this configuration. However, for large calculations, this runtime is the
most robust and reliable available.

If the communication network is Infiniband, you can configure GA to use the
Infiniband runtime instead by replacing `--with-mpi-pr` with
`--with-mpi --with-openib`. This runtime may be a little faster than
progress ranks. In this case the number of processors available to the GridPACK
application is the same as the number requested with `mpirun`. Also,
remember to set the `USE_PROGRESS_RANKS` variable to `FALSE` if
using this runtime.

## GridPACK

[Obtain the GridPACK release or development code](https://github.com/GridOPTICS/GridPACK/releases)
and put it in a convenient directory.  The
following commands can be used to configure and build GridPACK.

Make a directory in which to configure and build GridPACK (do not try to build
GridPACK in the `src` directory):

```
    cd path/to/gridpack/src
    mkdir build
    cd build
```

Use `cmake` to configure

```
    rm -rf CMake*

    cmake -Wdev \
        -D PETSC_DIR:STRING='/pic/projects/gridpack/software_new/petsc-3.7.6' \
        -D PETSC_ARCH:STRING='linux-openmpi-gnu-cxx-complex-opt' \
        -D GA_DIR:STRING='/pic/projects/gridpack/software_new/ga-5-6/build_pr' \
        -D USE_PROGRESS_RANKS:BOOL=TRUE \
        -D GA_EXTRA_LIBS='-lrt' \
        -D MPI_CXX_COMPILER:STRING='mpicxx' \
        -D MPI_C_COMPILER:STRING='mpicc' \
        -D MPIEXEC:STRING='mpiexec' \
        -D MPIEXEC_MAX_NUMPROCS:STRING="3" \
        -D BOOST_ROOT:STRING='/pic/projects/gridpack/software_new/boost_1_65_0' \
        -D CMAKE_INSTALL_PREFIX:PATH='/people/d3g293/gridpack/src/build_const/install' \
        -D CMAKE_BUILD_TYPE:STRING='RELWITHDEBINFO' \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        ..    
```

Alternatively, use the example configuration script (usually most up to date):

```
    sh ../example_configuration.sh constance
```

The argument to `PETS_ARCH` should be the same as the corresponding
string in the script to configure PETSc. This script assumes that ParMETIS is
part of the PETSc build. If ParMETIS was built separately, then you can include
the location of the ParMETIS libraries in the configuration script by adding the
line

```
    -D PARMETIS_DIR:STRING= '/path/to/parmetis/libraries' \
```

to the GridPACK configuration script.

This script assumes the GA libraries were built using the Infiniband runtime for
GA. The directories in this script need to be modified to reflect your system.
If the Infiniband runtime for GA is used, then `USE_PROGRESS_RANKS`
should be set to `FALSE`. The `MPIEXEC_MAX_NUMPROCS`,
which controls the maximum number of processors that are used when running the
GridPACK tests, should also be changed to 2.

After configuring GridPACK, build it with

```
    make
```

If compilation is successful, the [unit tests](BUILD_GRIDPACK.md) and/or
[example applications](BUILD_GRIDPACK.md) can be run.


## Issues

This build has been tested on the PNNL institutional cluster using the software
versions list at the top of this page. There are potential improvements that
could be made to this build:

* The GNU compilers can take advantage of the cluster hardware if certain options are used. They were not.
* The PETSc build should use the AMCL libraries for BLAS/LAPACK support rather than downloading CBlas.

