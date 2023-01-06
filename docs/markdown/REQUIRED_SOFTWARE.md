- [Building on Specific Platforms](#Building-on-Specific-Platforms)
- [CMake/CTest](#CMake/CTest)
- [MPI](#MPI)
- [Global Arrays](#Global-Arrays)
- [Boost](#Boost)
- [PETSc](#PETSc)
- [ParMETIS](#ParMETIS)
- [Doxygen](#Doxygen)
- [Linux Basics](#Linux-Basics)

## Building on Specific Platforms

Scripts or instructions for installing the libraries listed below on specific
platforms can be found by following the links below. We highly recommend that
you find an example that resembles the platform you are intending to use. Use
the scripts for that example as the basis for building GridPACK on your own
system. Most scripts require only minor modifications to get them working. These
usually involve changing file names so that they reflect your local directory
structure.

* ~~[Mac OS X (Yosemite) with MacPorts](DUMMY.md)~~ (outdated)
* [Mac OS X (High Sierra) with MacPorts](DUMMY.md)
* [Red Hat Enterprise Linux Workstation](DUMMY.md)
* [CentOS 6](DUMMY.md)
* [Ubuntu Linux](UBUNTU_LINUX.md)
* [Debian Linux](DUMMY.md)
* [PNNL PIC Cluster (Linux cluster with Infiniband)](DUMMY.md)
<!-- * [[Building_on_Windows | Native Windows Port]] -->

Notes on some additional builds can also be found in the docs/notes directory
under the top-level GridPACK source tree. These notes are for leadership class
parallel facilities that are much more complicated than most clusters. However,
they may have some useful ideas for handling problems if you run into
difficulties.

Many clusters use modules to install software such as the compilers, MPI
libraries, CMake, Git, etc. If you have such a system, you can install
compilers, MPI and CMake using a set of commands such as

```
   module purge
   module load gcc
   module load openmpi
   module cmake
```

These automatically install the libraries and modify environment variables so
that binaries and executables are in your path. Many systems may also have
modules for Boost and PETSc, but we urge users to be more cautious about using
these. These libraries are frequently not built with features that may be needed
by individual applications. For example, system builds of Boost often lack MPI
and PETSc is frequently missing the C++ interface. Users are probably better off
using their own library builds. This also guarantees that all libraries are
built with the same compiler and version of MPI. We suggest that you create a
library under your home directory that can be used to store all the libraries
used to build GridPACK and build and install all libraries (with the possible
exception of MPI) in this directory.

Some basic knowledge of Linux is also necessary for downloading libraries and
constructing build scripts. Information on using Linux is provided at the bottom
of this page. Additional information is available on the web.

## CMake/CTest

GridPACK uses the [cross-platform build system](http://www.cmake.org/)
cross-platform build system.  A
reasonably modern version should be used. Currently, we require version 2.8.8 or
above. You can check which version of CMake is on your machine by typing

```
  cmake -version
```

[CMake](http://www.cmake.org/) projects are designed to be built outside of the
source code location.  In the source directory of a GridPACK release
(GRIDPACK/src) create a subdirectory to use as the location of the build. (In
this documentation, GRIDPACK stands for the location of the top level GridPACK
directory.) The GRIDPACK/src directory should contain a file called
CMakeList.txt. Configure and build GridPACK in the subdirectory. The
`cmake` command should be pointed at the GRIDPACK/src directory that
contains the top-level CMakeList.txt file. If you are building in a directory
under GRIDPACK/src such as GRIDPACK/src/build, then your cmake command will have
the form

```
  cmake [ARGUMENTS] ..
```

The terminal `..` point the `cmake` command to the
directory immediately above the build directory. The `..` appear in
all example scripts listed in the builds on specific platforms. Make sure you
include them in your configure script!

Example configure scripts and other information on building CMake can be found
on the link below.

* [Redhat Linux Workstation](DUMMY.md)

For most systems, it is possible to install CMake using modules or an
installation capability such as `yum`.

## MPI

A working MPI implementation is required for building GridPACK. We commonly use
[OpenMPI](http://www.open-mpi.org/) and [MPICH](https://www.mpich.org). Recent
versions include OpenMPI 1.8.3 and MPICH 3.2. Other implementations have, such
as Intel MPI, also been used successfully.  Most MPI installations have compiler
wrappers, such as `mpicc` and `mpicxx` that combine the
compiler with all the directives needed to link to the MPI libraries. These can
be used directly in the GridPACK configuration.

Identify the compilers and `mpiexec` to the GridPACK configuration by
including `cmake` options like:

```
    -D MPI_CXX_COMPILER:STRING='mpicxx' \
    -D MPI_C_COMPILER:STRING='mpicc' \
    -D MPIEXEC:STRING='mpiexec' \
```

You can check to see if these wrappers are available on your system by typing

```
  which mpicc
```

If the wrapper is available, you should see a listing pointing to the location
of the wrapper. If you don't, you will need to modify your environment to
include it. Note that although `mpicc` and `mpicxx` are
fairly common names for the compiler wrappers, there is no standard and other
implementations of MPI may use something completely different. Check with your
system consultant for more information. Depending on the version of MPI you are
using, you may be able to find out more information by typing

```
  mpicxx -v
```

Other options may be needed by CMake to specify the MPI environment.  See the
documentation [http://www.cmake.org/cmake/help/v2.8.8/cmake.html#module:FindMPI
here]. An example configuration script for building MPI can be found on the link
below.

```
* [Redhat Linux Workstation](DUMMY.md)
```

On most systems it is possible to install MPI using modules or an installation
capability such as `yum`.

## Global Arrays

GridPACK depends heavily on [Global Arrays](https://github.com/GlobalArrays/ga).
The download page for recent releases can be found
[here](https://github.com/GlobalArrays/ga/releases).  The GA libraries used with
GridPACK must have the C++ interface enabled and the Fortran interface disabled.
More information on building GA can be found in the descriptions for building
GridPACK on individual platforms. The GridPACK configuration is not able to
identify additional required libraries if the Fortran interface is enabled or
independent BLAS/LAPACK libraries are used.  The following configuration options
should always be included when configuring GA on any platform

```
    --enable-cxx --without-blas --disable-f77
```

The `--without-blas` guarantees that GA does not try and build with the
BLAS libraries (which are downloaded and built with PETSc).

To configure GridPACK to recognize GA, specify the directory where Global Arrays
is installed and any extra libraries that are required:

```
    -D GA_DIR:PATH=/path/to/ga/install \
    -D GA_EXTRA_LIBS:STRING="..." \
    -D USE_PROGRESS_RANKS:BOOL=FALSE \
```

The `GA_EXTRA_LIBS` variable is used to include required libraries not
identified in the configuration. The `USE_PROGRESS_RANKS` variable
depends on the runtime used to build GA and should only be set to `TRUE`
if GA was configured using the `--with-mpi-pr` option.

We have used three different configuration of GA to build and run GridPACK. If
you are using GridPACK on a Linux cluster with an Infiniband interconnect, then
you can use the OpenIB runtime by including the `--with-openib` option
when configuring GA. This is the highest performing version of GA for clusters
with Infiniband, although for large calculations you can run into problems with
memory allocation. For any system with a working version of MPI, you can also
use the MPI two-sided runtime or the progress ranks runtime with GA. Use the
`--with-mpi-ts` or `--with-mpi-pr` options when configuring GA.
The two-sided runtime is the simplest runtime and is suitable for workstations
with a limited number of cores. This runtime provides reasonable performance on
a small number of cores but slows down considerably at larger core counts(our
experience is that you should limit this runtime to 8 or less processors). It is
not recommended for large-scale parallel computation.  The progress ranks
runtime is much higher performing and approaches the performance of the OpenIB
runtime. It is very reliable and runs on any platform that supports MPI.
However, it has one peculiarity in that it reserves one MPI process on each SMP
node to act as a communication manager. Thus, if you are running your
calculation on 2 nodes with 5 processes on each node, the GridPACK application
will only see 8 processes (4 on each node). To make sure that the GridPACK build
is aware of this, the `USE_PROGRESS_RANKS` parameter should be set to
`TRUE` when using the progress ranks build of GA.

Example scripts for configuring GA can be found in the links below.

* [Mac Yosemite ](DUMMY.md)
* [Mac High Sierra](DUMMY.md)
* [Redhat Linux Workstation](DUMMY.md)
* [CentOS 6](DUMMY.md)
* [Redhat Linux Cluster](DUMMY.md)

## Boost

The [Boost C++ Library](http://www.boost.org/) is used heavily throughout the
GridPACK framework, and a relatively recent version is required.  The GridPACK
configuration requires version 1.49 or later.  The
[Boost](http://www.boost.org/) installation must include
[Boost::MPI](http://www.boost.org/doc/libs/1_53_0/doc/html/mpi.html)
which must have been built with the same MPI compiler used for GridPACK.
Be aware that many Boost modules do not include MPI, so you may have to build
Boost on your own even if it is available as a module on your system.

To configure GridPACK to recognize Boost, one need only specify where
[Boost](http://www.boost.org/) is installed, like this

```
    -D BOOST_ROOT:STRING='/path/to/boost' \
```

Boost is tied quite closely to the latest features in C++ and problems can be
encountered if the version of Boost that you are using was released much later
than the compiler. Reverting to an earlier Boost version can sometimes eliminate
problems if you are having difficulties building it. The same is true for Boost
and CMake. If the CMake version was released earlier than the Boost version,
CMake may have problems identifying the libraries in Boost that it needs for
GridPACK. Again, going to an earlier version of Boost may fix these issues.

The Boost build can be tricky. Some clusters have Boost modules that can
potentially be used instead of building Boost on your own, but many modules are
not built with Boost::MPI. You will still need to specify the location of the
Boost directory, which can be found by using the command

```
    module show boost
```

This should tell you the location of `BOOST_ROOT`, which you can then use
in your GridPACK configuration script. If the GridPACK configuration does not
report that MPI was found, then you will need to get your system administrator
to rebuild Boost with the MPI libraries or build boost on your own. A successful
Boost configuration in GridPACK should report the results

```
    -- Checking Boost ...
    -- Boost version: 1.61.0
    -- Found the following Boost libraries:  
    --   mpi
    --   serialization
    --   random
    --   filesystem
    --   system
```

If you need to build Boost yourself, refer to the documentation on building
GridPACK on individual platforms for additional details on build Boost. If an
attempt to configure and build Boost fails, it usually is a good idea to fix the
build script and then remove the existing Boost directory and create a new one
by untarring the Boost tarball. Attempts to resume a failed Boost build after
fixing the build script are usually unsuccessful.

Example scripts for configuring and building Boost can be found on the links
below.

* [Mac Yosemite ](DUMMY.md)
* [Redhat Linux Workstation](DUMMY.md)
* [CentOS 6](DUMMY.md)
* [Redhat Linux Cluster](DUMMY.md)

## PETSc

GridPACK currently relies on the [Portable, Extensible Toolkit for Scientific
Computation (PETSc)][http://www.mcs.anl.gov/petsc/index.html) for parallel
linear algebra, and linear and nonlinear system solvers. The PETSc interface
tends to change a bit as new releases come out, requiring adjustments in any
applications that use it. We have currently used PETSc versions 3.4-3.8 with
GridPACK.

PETSc is a complicated package with numerous options.  PETSc needs to be built
with MPI enabled and using the same MPI implementation used for GridPACK and the
other libraries such as Boost and GA.  It also needs to use C++ as the base
language. Originally, GridPACK could only use PETSc if it was configured for
complex support. The current GridPACK release can use either complex or real
builds. However, most applications in GridPACK use complex matrices, so it is
still preferable to configure PETSc to use complex variables. Refer to the
[PETSc installation documentation](http://www.mcs.anl.gov/petsc/documentation/installation.html)
for additional information on how to configure PETSc.

Configuring and building PETSc is done in the top level PETSc directory. One of
the configuration variables that needs to be set when configuring and building
PETSc is PETSC_ARCH. In the example below, PETSC_ARCH was set to
`'arch-Darwin-cxx-opt'`. After the build is complete, there will be a
directory beneath the top level directory with whatever name was assigned to
PETSC_ARCH. This directory contains the include and lib directories for the
PETSc libraries.

The GridPACK configuration must know where
[PETSc](http://www.mcs.anl.gov/petsc/index.html) is installed.  This is specified
by two options as shown below. 

    -D PETSC_DIR:STRING='/Users/d3g096/ProjectStuff/petsc-3.4.0' \
    -D PETSC_ARCH:STRING='arch-darwin-cxx-opt' \

Currently, the configuration will recognize and adjust the GridPACK build if the
[PETSc](http://www.mcs.anl.gov/petsc/index.html) build includes
[ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview),
[SuperLU_DIST](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) and/or
[MUMPS](http://mumps.enseeiht.fr/).  Many of the example GridPACK applications
expect a parallel direct linear solver to be built into PETSc.  This is
satisfied by including
[SuperLU_DIST](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/)
or [MUMPS](http://mumps.enseeiht.fr/) in the PETSc build.

Examples of complete scripts for configuring PETSc can be found on the links
below.

* [Mac Yosemite ](DUMMY.md)
* [Mac High Sierra](DUMMY.md)
* [Redhat Linux Workstation](DUMMY.md)
* [CentOS 6](DUMMY.md)
* [Redhat Linux Cluster](DUMMY.md)

## ParMETIS

GridPACK uses
[ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
to (re)distribute a power grid network over several processors.  It
needs to be built with the same MPI configuration as
[Boost](http://www.boost.org/)
and [PETSc](http://www.mcs.anl.gov/petsc/index.html). The GridPACK
configuration script will find
[ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
automatically if it has been included in the
[PETSc](http://www.mcs.anl.gov/petsc/index.html) build. Otherwise, the GridPACK
configuration just needs to know where
[ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview) was
installed, which is specified by

```
    -D PARMETIS_DIR:STRING='/pic/projects/gridpack/software' \
```

GridPACK requires version ParMETIS version 4.0.  Older versions will not work.
On most systems, it is straightforward to download and build ParMETIS as part of
the PETSc build. We highly recommend that you do this to access ParMETIS.

An example configuration script for ParMETIS can be found on the link below.

*[Redhat Linux Workstation](DUMMY.md)

## Doxygen

GridPACK uses [Doxygen](http://www.doxygen.org/) to help document code. It's use
is optional. [Doxygen](http://www.doxygen.org/) documentation can optionally be
prepared during the build process.  This is enabled if
[Doxygen](http://www.doxygen.org/)
is found.  [Graphviz](http://www.graphviz.org/) is necessary for full
documentation features.

## Linux Basics

Some basic knowledge of Linux is necessary in order to build GridPACK.
Familiarity with a Linux editor such as VIM or EMACS is required. Extensive
documentation is readily available for both these editors, both online and in
books. In addition, you will need to download tarballs (.tar or .tar.gz files)
of the Boost, PETSc and GA libraries, uncompress them, and then configure and
build the libraries.

The most direct way to download tar files on a Linux machine is to use a web
browser such as Firefox to download files to a local download directory. These
can then be moved to another location such are your GridPACK software directory.

If, for some reason, you do not have a web browser on your machine, you can use
the `wget` command. This has the form

```
    wget URL
```

where URL is the location of a tar file on the web. The difficult part of using
this command is to identify the URL of the download file. You will still need to
have a web browser working on some other platform that you can use to navigate
to the download page. If you go to the download page for these libraries and
mouse over one of the download links for a tar file (these have the extension
`.tar.gz`), you will see the URL of the download in the lower left hand
corner. This seems to work for several different web browsers. Using this URL in
the wget command will download the tar file into your current directory. This is
the simplest method for downloading tar files.

For Windows Explorer, if you mouse over the link and right-click, a menu pops
up. Select the "Copy shortcut" option. If you are accessing the Linux machine
via a Putty session on your Windows box, then you can use a right click after
typing `wget` in the command prompt to paste the URL directly into your
Linux window.

If you usually download files using a Windows machine, you can download a
library tar file using Windows and then copy it to your Linux platform using the
[https://winscp.net/eng/download.php WinSCP] utility. This will allow you to
transfer your files from Windows to Linux in a straightforward way. You may also
be able to download directly to a Linux platform by bringing up a browser such
as Firefox from the Linux command prompt, going to the appropriate download site
and downloading directly to your Linux directory.

Once you have a tarball downloaded to your software directory, the next step is
to uncompress the file into its own directory. For example, if you have
downloaded the Boost tarball for version 1.65.0, you would see the following
file in your directory

```
    boost_1_65_0.tar.gz
```

The .tar extension means that all the files in the boost directory have been
concatenated into a single file using the Linux tar command. The .gz extension
means that the tarball has been further compressed using the gzip command. You
can uncompress the file and untar it using the single command

```
    tar xvf boost_1_65_0.tar.gz
```

This will produce a directory

```
    boost_1_65_0
```

in the same directory that the tarball is located in. This directory will
contain all the individual files and subdirectories in the Boost library. Other
libraries have a similar structure.

Continuing with the Boost example, the next step is to cd into the Boost
directory and create a script for building Boost. This would consist of creating
a file at the top of the boost_1_65_0 directory with a name such as
`build.sh` that contains the commands for configuring Boost. On a
Redhat Linux cluster using the GNU compilers, you would use these lines

```
    echo "using mpi ;" > ~/user-config.jam
    sh ./bootstrap.sh \
        --prefix="/my_home_directory/software/boost_1_65_0" \
        --without-icu \
        --with-toolset=gcc \
        --without-libraries=python
    ./b2 -a -d+2 link=static stage
    ./b2 -a -d+2 link=static install
    rm ~/user-config.jam
```

Note that the argument to `--prefix` is the path to the Boost
directory that you are currently in. Once these lines have been copied into the
`build.sh` file, the file needs to be made executable by changing its
permissions with the command

```
    chmod +x build.sh
```

The script can then be run (in the boost_1_65_0 directory) by typing

```
    ./build.sh
```

This will configure and build Boost. It is not strictly necessary to put these
commands in a script, they will also work by just typing them into the the Linux
command line. However, for such a long set of commands, it is obviously more
desirable to avoid a mistake by using the script.

A final note on Linux concerns which shell is being used in your window session.
For most users this is either the C-shell or Bourne-shell. If you are using
C-shell, then environment variables are set using a command of the form

```
    setenv FC gfortran
```

where `FC` is the environment variable and `gfortran` is the
value. For the Bourne-shell, the command is different

```
    export FC=gfortran
```

This is the only distinction between the two shells that is needed when building
GridPACK and its associated libraries. All other scripts can be used in either
shell without modification.

