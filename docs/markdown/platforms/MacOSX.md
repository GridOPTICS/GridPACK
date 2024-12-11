Introduction
============

This details the building of GridPACK on a Macbook Pro running Mac OS X
Sonoma (14.7.1). The build was aided by the use of
[MacPorts](https://www.macports.org/) to install compilers and MPI.
Super user permissions are required to install and manage
[MacPorts](https://www.macports.org/).

These instructions use a specific combination of compilers, MPI
implementation, and [MacPorts](https://www.macports.org/) packages.
Other combinations can be probably used, the one here happened to work
for the author. If you choose some other combination, you\'re on your
own.

The author uses the C-shell. These instructions will need to be slightly
modified to adapt to another shell, like `bash`.

MacPorts
========

[MacPorts](https://www.macports.org/) is free, open source system
whereby various software packages can be installed and maintained on a
Mac OS X system. [Homebrew](https://brew.sh/) provides a similar
capability, but the author has no experience with it.

The [MacPorts](https://www.macports.org/) project provides very good
[installation instructions](https://www.macports.org/install.php). They
will not be repeated here. It is assumed here that
[MacPorts](https://www.macports.org/) was installed in `/opt/local`,
which is the default.

After [MacPorts](https://www.macports.org/) installation,
`/opt/local/bin` needs to be in shell search path. The easiest way to do
this is to add `/opt/local/bin` to the top of `/etc/paths` (superuser
privileges required). Afterward, a reboot or starting a new shell should
make the `port` command availble on the command line.

Packages from MacPorts
----------------------

Install the CLang 16.0 compiler set and mpich MPI implementation.

    sudo port install mpich mpich-clang16
    sudo port select --set clang mp-clang-16
    sudo port select --set mpi mpich-clang16-fortran
    sudo port select gcc mp-gcc14

The latter commands make the installed compilers, and MPI wrappers the
default. That way the compiler commands, `clang` and `clang++`, and MPI
compiler wrappers, like `mpicc`, will be available at the command line.

Python 3 is necessary, so a specific version was chosen. Whatever
version should be adequate.

    sudo port install python312
    sudo port select --set python python312
    sudo port select --set python3 python312

Install CMake and other necessities:

    sudo port install cmake
    sudo port install pkgconfig


Build GridPACK with Included Scripts
====================================

The `install_gridpack_deps.sh` and `install_gridpack.sh` scripts
included in the GridPACK repository can usually simplify the GridPACK
build process.  On this example system, using the MacPorts packages as
described above, GridPACK was successfully built using those scripts
as follows (must be executed within a `bash` shell:

    bash
    git clone https://github.com/GridOPTICS/GridPACK.git
    cd GridPACK/
    export LD_LIBRARY_PATH=/opt/local/lib/gcc14
    source install_gridpack_deps.sh
    export CC=/opt/local/bin/clang
    export CXX=/opt/local/bin/clang++
    source install_gridpack.sh
    cd src/build
    ctest
    cd applications/dynamic_simulation_full_y
    ../../../install/bin/dsf2.py input_9b3g.xml

Build from Source
=================

Prerequisites
-------------

[PETSc](http://www.mcs.anl.gov/petsc/) and [Global
Arrays](http://www.emsl.pnl.gov/docs/global/) need to be built from
source. You need to choose a path in which packages will be installed.
Here the author decided to use `$HOME/Projects/GridPACK/install` as an
installation directory. Below, it is assumed that `$prefix` is set

    set prefix = $HOME/Projects/GridPACK/install

In order to build the GridPACK Python interface, all prerequisites must
be built as shared libraries.

I tried to install Boost and PETSc via MacPorts. Boost, with proper
variants, will install and is usable to GridPACK. However, PETSc 3.22
could not be built with SuperLU or Suitesparse. This tends to break some
solver unit tests. So, build from source.

### Global Arrays ###

GA version 5.8.2 source code was downloaded from the [Github
site](https://github.com/GlobalArrays/ga/releases/download/v5.8.2/ga-5.8.2.tar.gz).
Avoid cloning the Github repository, because preparing to build the code
requires a full installation of GNU autotools (and that has been
problematic for me).

GA was configured and built as follows

    cd ga-5.8.2
    mkdir build
    cd build
    unsetenv CC
    unsetenv CXX
    ../configure --with-mpi-ts --disable-f77 \
        --without-blas --without-lapack --without-scalapack \
        --enable-cxx --enable-i4 \
        --prefix=$prefix \
        --enable-shared=yes --enable-static=no \
        CFLAGS=-Wno-implicit-function-declaration
    make -j
    make -k check
    make install

### Boost ###

Boost version 1.81.0 was downloaded from
[here](https://archives.boost.io/release/1.81.0/source/).

In the past, it\'s tricky to get Boost build to use the MacPorts CLang,
but it seems to work here. The build uses \"clang++\" on PATH rather
than a full path.

    ./bootstrap.sh --with-toolset=clang --prefix=$prefix \
        --with-libraries=mpi,serialization,random,filesystem,system
    echo 'using mpi : mpicxx ; ' >> project-config.jam
    ./b2 -a -d+2 toolset=clang link=shared stage
    ./b2 -a -d+2 toolset=clang link=shared install

### PETSc ###

On modern Mac OS X, only PETSc version 3.20.1 and later can be used.
Older versions will not build. This has been noted by others.

    wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.20.6.tar.gz
    tar xvzf petsc-3.20.6.tar.gz
    cd petsc-3.20.6
    unsetenv PETSC_ARCH
    python ./configure \
        PETSC_ARCH=mpich-clang-real-shared \
        --prefix=$prefix \
        --with-mpi=1 \
        --with-gnu-compilers=0 \
        --with-cc=/opt/local/bin/mpicc \
        --with-fc=/opt/local/bin/mpif90 \
        --with-cxx=/opt/local/bin/mpicxx \
        --with-clanguage=c \
        --with-fortran-bindings=0 \
        --with-scalar-type=real \
        --with-precision=double \
        --download-cmake=0 \
        --download-suitesparse=1 \
        --download-superlu_dist=1 \
        --download-superlu=1 \
        --download-parmetis=1 \
        --download-metis=1 \
        --download-f2cblaslapack=1 \
        --download-mumps=1 \
        --download-scalapack=1 \
        --with-shared-libraries=1 \
        --with-static-libraries=0 \
        --with-x=0 \
        --with-valgrind=0 \
        --with-mpiexec=/opt/local/bin/mpiexec \
        --with-debugging=0
    make PETSC_DIR=/Users/d3g096/Projects/GridPACK/src/petsc-3.20.6 \
        PETSC_ARCH=mpich-clang-real-shared all
    make PETSC_DIR=/Users/d3g096/Projects/GridPACK/src/petsc-3.20.6 \
        PETSC_ARCH=mpich-clang-real-shared install

Building GridPACK
-----------------

CMake will work very hard to choose the wrong compiler. Make sure it
uses the correct compiler by setting environment variables.

Configure, build, and test GridPACK with

    unsetenv DYLD_LIBRARY_PATH
    setenv CC /opt/local/bin/clang
    setenv CXX /opt/local/bin/clang++
    rm -rf CMake*
    cmake \
        -D BUILD_GA:BOOL=NO \
        -D CMAKE_BUILD_TYPE:STRING=Release \
        -D GA_DIR:PATH=$prefix \
        -D BOOST_ROOT:STRING=$prefix \
        -D Boost_NO_BOOST_CMAKE:BOOL=TRUE \
        -D PETSC_DIR:PATH=$prefix \
        -D MPI_CXX_COMPILER:STRING='/opt/local/bin/mpicxx' \
        -D MPI_C_COMPILER:STRING='/opt/local/bin/mpicc' \
        -D MPIEXEC:STRING='/opt/local/bin/mpiexec' \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=120 \
        -D USE_GLPK:BOOL=OFF \
        -D BUILD_SHARED_LIBS:BOOL=YES \
        -D CMAKE_INSTALL_PREFIX:PATH=$prefix \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=Yes \
        ..
    make -j 4
    ctest
    make install

GridPACK Python Interface
-------------------------

An additional MacPorts package is needed for the Python interface

    sudo port install py312-mpi4py +clang16 +mpich
    sudo port install py-pip
    sudo port select --set pip pip312
    sudo port select --set pip3 pip312

Alternatively, a Python virtual environment can be used to install
`mpi4py` if desired.

The Python version chosen here is new enough that it does not like to
run `setup.py` directly, as directed in README.md. So, use `pip` to
install the module:

    setenv GRIDPACK_DIR $prefix
    cd .../GridPACK/python
    pip install --no-deps --upgrade --prefix=$prefix .

This will install the Python module in the same place as the rest of
GridPACK (`$GRIDPACK_DIR`). To use the module, `PYTHONPATH` must be set.
In this case, it had to be

    setenv PYTHONPATH $prefix/lib/python3.12/site-packages/

An easy check to see if the module works is

    python -c 'import gridpack'

The first attempt resulted in

    Traceback (most recent call last):
      File "<string>", line 1, in <module>
    ImportError: dlopen(/Users/d3g096/Projects/GridPACK/install/lib/python3.12/site-packages/gridpack.cpython-312-darwin.so, 0x0002): Library not loaded: @rpath/libgfortran.5.dylib
      Referenced from: <B4CFCC4B-7ACE-3880-BA09-FFC95B7F583E> /Users/d3g096/Projects/GridPACK/install/lib/python3.12/site-packages/gridpack.cpython-312-darwin.so
      Reason: tried: '/Users/d3g096/Projects/GridPACK/install/lib/libgfortran.5.dylib' (no such file), '/System/Volumes/Preboot/Cryptexes/OS/Users/d3g096/Projects/GridPACK/install/lib/libgfortran.5.dylib' (no such file), '/Users/d3g096/Projects/GridPACK/install/lib/libgfortran.5.dylib' (no such file), '/System/Volumes/Preboot/Cryptexes/OS/Users/d3g096/Projects/GridPACK/install/lib/libgfortran.5.dylib' (no such file)

For some reason, the Fortran shared library could not be found, probably
something in the PETSc build, so the `DYLD_LIBRARY_PATH` environment
variable needs to be set

    setenv DYLD_LIBRARY_PATH /opt/local/lib/gcc14

Some complete Python script tests can be found in
`.../python/src/example`. These exercise the HADREC and dynamic
simulation applications. Once `PYTHONPATH` and `DYLD_LIBRARY_PATH` are
set, these should complete without error

    cd src/example
    python 39bus_test_example.py
    python 39bus_test_example_dsf.py
    python 39bus_scatterload_steptest_new_itr.py
    mpiexec -np 2 python 39bus_scatterload_steptest_new_itr.py
    python 39bus_scatterload_steptest_new_itr_dsf.py
    mpiexec -np 2 python 39bus_scatterload_steptest_new_itr_dsf.py
    python 39bus_scatterload_steptest_new_itr_compensateY.py
    python 39bus_test_pfdata.py

