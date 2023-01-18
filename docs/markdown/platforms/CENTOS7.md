This documents a GridPACK build on [CentOS](https://www.centos.org/ CentOS) 7.5.1804
installed on a [VirtualBox](https://www.virtualbox.org/) virtual machine (VM).
Virtual Box is not required and these instructions should work for any Linux
system using the CentOS operating sytem.  They also work for
[Red Hat Enterprise Linux](https://www.redhat.com/en/technologies/linux-platforms/enterprise-linux)
(RHEL) 7.

== Requisite System Software

In this build, some [prerequisite software](../required/BASIC_LINUX.md) was installed
as system software packages. If you are using a workstation with a standard
configuration, there is a good chance that this software is already installed
but if you are setting up a VM then you may need to add this software on your
own.  *Super user priviledges are required to install system
software.* If you do not have such priviledges, you will need a system
administrator to perform some of these steps. Super user priviledges are not
required to build and use GridPACK, only to install system software.

Most required system software is available in the default CentOS/RHEL
repostitories. However, at least one package from the
[Extra Packages for Enterprise Linux](https://fedoraproject.org/wiki/EPEL)
(EPEL) repository is needed.  Add access to the EPEL repository with

```
    sudo yum install epel-release
```

Install a C++ compiler, [CMake](https://cmake.org/), and
[Git](https://git-scm.com/) with

```
    sudo yum install git cmake cmake3 gcc-c++
```

`cmake` is required to build some of the packages. `cmake3` is required for
GridPACK.

=== OpenMPI

Install OpenMPI (version 1.10) using

```
    sudo yum install openmpi openmpi-devel
```

This will install the `module` command. Before continuing, it will be
necessary to either open another terminal and proceed there, or log out and log
in. This makes the `module` command available for use. OpenMPI is
installed as a "module". In order to use `mpiexec` and compiler
wrappers, the correct module must be loaded using

```
    module load mpi/openmpi-x86_64
```

This makes the MPI compiler wrappers and `mpiexec` available to the
command line. *The OpenMPI module must be loaded each
time a terminal session is started.* If you are using OpenMPI across multiple
projects, you may want to put the module load command in `.cshrc` or `.bashrc` file
so that it is automatically loaded whenever you create a new terminal.

=== Boost

Install [Boost](https://www.boost.org/) (version 1.53) with

```
    sudo yum install boost boost-devel boost-openmpi boost-openmpi-devel
```

This will install the needed Boost libraries and headers.

== Build Packages from Source

In this build, all packages built from source (and GridPACK) were installed in a
common directory, `$HOME/gridpack`.
For convenience, a shell variable was set to this directory using, for the
C-shell

```
    set PREFIX = $HOME/gridpack<
```

or in the Bourne shell

```
    export PREFIX=$HOME/gridpack
```

The common directory appears as `$PREFIX` below. 

Building GridPACK on CentOS/RHEL 7 requires code be built with pthreads
throughout. Before continuing, it will be necessary to set environment variables
to make sure proper compilers and their pthread options are used. If using the
C-shell:

```
    setenv CC gcc
    setenv CFLAGS "-pthread"
    setenv CXX g++
    setenv CXXFLAGS "-pthread"
    setenv FC gfortran
    setenv FFLAGS "-pthread"
```

or in a Bourne shell

```
    export CC=gcc
    export CFLAGS="-pthread"
    export CXX=g++
    export CXXFLAGS="-pthread"
    export FC=gfortran
    export FFLAGS="-pthread"
```

=== PETSc, version 3.16.3

A PETSc package is available through `yum` (EPEL), but its headers
and libraries are installed in locations that the GridPACK configuration does
not understand. It's best to just build PETSc from source.

Get the PETSc 3.16.3 source from the
[PETSc download page](https://www.mcs.anl.gov/petsc/download/index.html) and
unpack it. In this case, it was unpacked in `$HOME/gridpack/src/petsc-3.16.3`. In the
unpacked directory, configure PETSc as follows:

    ./configure \
       PETSC_ARCH=arch-linux2-complex-opt \
       --with-fortran-kernels=1 \
       --download-superlu_dist \
       --download-superlu \
       --download-parmetis \
       --download-metis \
       --download-suitesparse \
       --download-fblaslapack \
       --with-clanguage=c++ \
       --with-shared-libraries=0 \
       --with-scalar-type=complex \
       --with-x=0 \
       --with-mpiexec=mpiexec \
       --prefix="$PREFIX/petsc-3.16.3" \
       --with-debugging=0 CFLAGS=-pthread CXXFLAGS=-pthread FFLAGS=-pthread

This will build some required packages and configure PETSc for compilation. It
will complain about how the `CC`, etc. environment variables were
set, but it ignores those and does the right thing.  Build the PETSc libraries
using

```
    make PETSC_DIR=$HOME/gridpack/src/petsc-3.16.3 PETSC_ARCH=arch-linux2-complex-opt all
```

Install PETSc using

```
    make PETSC_DIR=$HOME/gridpack/src/petsc-3.16.3 PETSC_ARCH=arch-linux2-complex-opt install
```

The installation can be tested with

```
    make PETSC_DIR=$PREFIX/petsc-3.16.3 PETSC_ARCH="" test
```

=== Global Arrays, 5.8

The EPEL repository also has a [Global Arrays](http://hpc.pnl.gov/globalarrays/)
package that can be installed with `yum`. However, at the time of
writing, this package was not compatible with the GridPACK configuration. So,
building GA from soure is necessary. 

Get the [Global Arrays](http://hpc.pnl.gov/globalarrays/) source from the
[release page](https://github.com/GlobalArrays/ga/releases/tag/v5.8.2).  Unpack
in a convenient location. In the unpacked directory, configure, build, and
install using

```
    ./configure --with-mpi-ts --disable-f77 --without-blas \
        --enable-cxx --enable-i4 --prefix="$PREFIX"
    make
    make install
```

== GridPACK

[Obtain the GridPACK release or development code](https://www.gridpack.org/wiki/index.php/Download_GridPACK)
and put it in a convenient directory, like `$HOME/gridpack/src`. 

Then, change into the top-most GridPACK source directory, make a build directory
e.g.:

```
    cd $HOME/gridpack/src/GridPACK
    mkdir src/build
    cd src/build
```

Configure and build GridPACK with

```
    cmake3 -Wdev --debug-trycompile \
        -D GA_DIR:STRING="$PREFIX" \
        -D USE_PROGRESS_RANKS:BOOL=OFF \
        -D BOOST_ROOT:STRING="/usr" \
        -D PETSC_DIR:STRING="$PREFIX/petsc-3.7.7" \
        -D MPI_CXX_COMPILER:STRING="/usr/lib64/openmpi/bin/mpicxx" \
        -D MPI_C_COMPILER:STRING="/usr/lib64/openmpi/bin/mpicc" \
        -D MPIEXEC:STRING="/usr/lib64/openmpi/bin/mpiexec" \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=400 \
        -D CMAKE_INSTALL_PREFIX:PATH="$PREFIX" \
        -D BUILD_SHARED_LIBS:BOOL=NO \
        ..
    make
```

If compilation is successful, the
[unit tests and/or example applications](../required/GRIDPACK.md#running-tests)
can be run.
