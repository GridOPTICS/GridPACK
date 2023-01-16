## Modules

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
reasonably modern version should be used. Most Linux platforms come with a
version of CMake installed, but the version may not be high enough.
Currently, we require a version of 3.5.0 or
above. You can test for CMake on your system by typing

```
  which cmake
```

You should see something like

```
  /usr/bin/cmake
```

If you get a message saying that command is not found, then CMake is not on your
system and you will need to build it.
You can check which version of CMake is on your machine by typing

```
  cmake -version
```

Again, if the version is below 3.5.0, you may need to build your own version of
CMake.

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

For most systems, it is possible to install CMake using modules or an
installation capability such as `yum`.  In the case that this is not possible
or your version is too old, you can build it using the
commands below.

You will need
to start by downloading the CMake tar file (e.g. cmake-3.24.3.tar.gz) from the
[download site](https://cmake.org/download/). Information on downloading tar
files can be found [here](#linux-basics) (at the bottom of this page).

```
  tar xvzf cmake-3.24.3.tar.gz
  cd cmake-3.24.3
  ./bootstrap --prefix=$PREFIX
  make
  make test
  make install
```

The `$PREFIX` represents the location that you want to put the CMake
executables.

## MPI

A working MPI implementation is required for building GridPACK. We commonly use
[OpenMPI](http://www.open-mpi.org/) and [MPICH](https://www.mpich.org).
Other implementations, such as Intel MPI, have
also been used successfully.  Most MPI installations have compiler
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
documentation
[here](http://www.cmake.org/cmake/help/v2.8.8/cmake.html#module:FindMPI).

On most systems it is possible to install MPI using modules or an installation
capability such as `yum`.

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

