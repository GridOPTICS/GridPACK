## Obtaining GridPACK Source

GridPACK source code can be obtained in two ways. Users interested in obtaining
the latest improvements to GridPACK can download the source code directly by
cloning the [GridPACK Github repository](https://github.com/GridOPTICS/GridPACK).
The GridPACK build system uses some third-party code that is included as a
submodule. After cloning GridPACK, it is absolutely *imperative* to do the
following in the top-most directory of your clone to get the submodule code:

```
    git submodule update --init
```

The code in the repository is undergoing constant revisions and may not be in a
stable state. For most users, it is better to download the latest release from
the [GridPACK release page](https://github.com/GridOPTICS/GridPACK/releases).
This is in a compressed `tar` archive format file.   Unpack the
archive in a convenient location using something like  

```
    tar xvf gridpack-X.X.tar.gz
```

For further information on using `tar` archives, there are many
[resources online](https://www.howtogeek.com/248780/how-to-compress-and-extract-files-using-the-tar-command-on-linux/)
to help.

You do not need to run the git submodule command if you download one of the
release tarballs. The tarballs already have the third party modules included.

## Configuring GridPACK

Configuration is the most complicated part of the process of actually building
GridPACK. Before attempting to configure and build GridPACK, make sure that all
the [prerequisite software](PRERQUISITE.md) needed by GridPACK is in place.
The information in this section is a general overview for configuring
GridPACK and adjustments may be required to get the build to work on particular
platforms. [CMake](https://cmake.org/download) is used to configure GridPACK for
building. CMake is usually available on most Linux systems but older versions of
Linux may have a version of Cmake that is too old to build GridPACK. GridPACK
requires version 3.5.0 or newer. If the installed version of CMake is too old,
users will need to download CMake and build it themselves. The version of CMake
can be found by typing

```
  cmake -version
```

The configure process insures that required software is available and usable.
CMake expects to configure GridPACK in a directory other than the one containing
the source code.  Typically, one makes an empty directory, usually called
something like `build` and then cd's into this directory. The
commands are

```
  mkdir build
  cd build
```

Once in the build directory, execute the command 

```
  cmake [options] gridpack/source/directory
```

where `options` are used to locate required software
and set compiler options. The shell script
`example_configuration.sh` shows some examples of configuration
options for a few systems. If you don't get the configure right the first time,
then you should make sure that you get rid of all the files that CMake created
when you tried configuring previously. This can be done by typing

```
  rm -rf CMake*
```

in your build directory. This will remove all CMake-related configuration files
from your build directory so that the new build is not corrupted by the previous
build. If you are using a script to configure GridPACK, then you should include
this line at the start of the script.

To guarantee that CMake finds the correct C and C++ compilers, you should define
the environment variables CC and CXX (if you are trying to build the Fortran
interface, you should define FC and F77 as well). In addition, it may also be
necessary to define `CFLAGS = "-pthread"` etc. depending on how some of the other
libraries were built. Using C-shell, the environment variables are

```
  setenv CC gcc
  setenv CFLAGS "-pthread"
  setenv CXX g++
  setenv CXXFLAGS "-pthread"
  setenv FC gfortran
  setenv FCFLAGS "-pthread"
  setenv F77 gfortran
  setenv F77FLAGS "-pthread"
```

A complete configuration line for GridPACK is

```
rm -rf CMake*

cmake -Wdev \
      -D PETSC_DIR:STRING='/projects/gridpack/software/petsc-3.16.3' \
      -D PETSC_ARCH:STRING='linux-openmpi-gnu-cxx-complex-opt' \
      -D GA_DIR:STRING='/projects/gridpack/software/ga-5.7/build_pr' \
      -D USE_PROGRESS_RANKS:BOOL=TRUE \
      -D MPI_CXX_COMPILER:STRING='mpicxx' \
      -D MPI_C_COMPILER:STRING='mpicc' \
      -D MPIEXEC:STRING='mpiexec' \
      -D MPIEXEC_MAX_NUMPROCS:STRING="3" \
      -D BOOST_ROOT:STRING='/projects/gridpack/software/boost_1_65_0' \
      -D CMAKE_INSTALL_PREFIX:PATH='/people/d3g293/gridpack/src/build_pr/install' \
      -D CMAKE_BUILD_TYPE:STRING='DEBUG' \
      -D CMAKE_VERBOSE_MAKEFILE:STRING=TRUE \
      ..
```

This example assumes that the build directory is immediately below the
`GRIDPACK/src` directory so that the `..` at the end of this
command are pointing to `GRIDPACK/src`. The directories that are arguments to options like
`BOOST_ROOT` and `PETSC_DIR` etc. should be replaced with
locations of these libraries on your system. The `PETSC_ARCH` variable is
also likely to depend on the particular build and will need to be modified by
the user to reflect their system. This example assumes that PETSc was not built
using the `--prefix` option. If the `--prefix` option was used, then the
`PETSC_ARCH` variable does not need to be specified. The `rm -rf CMake*` command
before calling `cmake` clears out the configuration files from previous builds.
It is generally a good idea to get rid of these files, since they may cause
problems is there was an error on the previous attempt to build GridPACK.

If ParMETIS has been
downloaded and built as part of PETSc, it is unnecessary to specify the location
of the ParMETIS installation and the `PARMETIS_DIR` option does not need
to be included. If it ParMETIS has been installed separately, then the location
of the ParMETIS and METIS libraries is specified by `PARMETIS_DIR`.

The Boost install directory is specified with the `BOOST_ROOT` option.
Typically, this is all that needs to be done for Boost, but be aware that if
Boost was not built with the mpi, serialization, random, system and filesystem
libraries, the configuration may fail on Boost. The mpi library is frequently
not included in system library builds and modules, so it may be necessary for
users to build their own version of Boost, even if a library is available on the
system.

The installation directory for Global Arrays can be specified using the `GA_DIR`
option. Extra libraries and compiler flags that may be needed by GA can be
specified using the `GA_EXTRA_LIBS` option. For MPI-based runtimes, this
option is usually not needed. The example above assumes that GA was built with
the high performance progress ranks runtime
(`--with-mpi-pr`). In this case, the `USE_PROGRESS_RANKS` option must be
included and set to true. If the two-sided runtime in GA is used, then the
`USE_PROGRESS_RANKS` option can be left out or set to false.

The MPI compilers and execute commands are specified in the <tt>cmake</tt>
command with the options

```
    -D MPI_CXX_COMPILER:STRING='mpicxx' \
    -D MPI_C_COMPILER:STRING='mpicc' \
    -D MPIEXEC:STRING='mpiexec' \
```

These options set the MPI wrappers for both the C and C++ compilers.

Some GridPACK-specific options are available for controlling the configuration:

More information on the CMake options is
available in the GridPACK [overview document](docs/user_manual/GridPACK.pdf)
(see the section on configuring and building GridPACK). Additional
examples can also be found in the pages describing builds on different
platforms.

## Building

Once configured, GridPACK is built with

```
  make
```

which will take some time. (The make command should be available by default on
most Linux platforms.)  If building on a multi-core system, building can go
faster if multiple cores are use, e.g.,

```
  make -j 8
```

will use 8 simultaneous processes to build GridPACK. It is possible that you may
run into problems building with multiple cores since some dependencies may get
out of order. If the build stops, then try continuing using only 1 core.

## Running Tests

After a successful build, GridPACK unit tests can be run with

```
  make test
```

or
 
```
  ctest
```

at the top of the GridPACK build directory. This will produce a list of tests
and whether they passed or failed. While every effort is taken to make sure
these tests pass, several things can lead to test failure, not just errors in
GridPACK code, including OS, compilers, and requisite software versions. We
cannot test all permutations.  Typical output from `ctest` might look
like this (abridged example): 

```
    Running tests...
    /usr/bin/ctest3 --force-new-ctest-process 
    Test project
/home/d3g096/Projects/GridPakLDRD/GridPACK.github/master/src/build
          Start  1: greetings_serial
     1/87 Test  #1: greetings_serial .....................   Passed    0.27 sec
          Start  2: greetings_parallel
     2/87 Test  #2: greetings_parallel ...................   Passed    0.32 sec
    [...]
          Start 68: state_estimation_serial
    68/87 Test #68: state_estimation_serial ..............   Passed    0.30 sec
          Start 69: state_estimation_parallel
    69/87 Test #69: state_estimation_parallel ............***Failed    0.33 sec
          Start 70: kalman_ds_serial
    70/87 Test #70: kalman_ds_serial .....................   Passed    2.88 sec
          Start 71: kalman_ds_parallel
    71/87 Test #71: kalman_ds_parallel ...................***Failed    0.45 sec
          Start 72: hello_world_serial
    72/87 Test #72: hello_world_serial ...................   Passed    0.29 sec
          Start 73: hello_world_parallel
    73/87 Test #73: hello_world_parallel .................   Passed    0.34 sec
    [...]
    
    98% tests passed, 2 tests failed out of 87
    
    Total Test time (real) =  67.29 sec
    
    The following tests FAILED:
             69 - state_estimation_parallel (Failed)
             71 - kalman_ds_parallel (Failed)
    Errors while running CTest
```

In most cases, there is a `serial` and `parallel` versions
of each test.  Note that, in this example, two tests were reported as failed and
the reason was very helpfully reported as `Failed`.  If the test is
marked as `Failed` without any further information, it could mean the
test crashed and did not run to completion or it did run to completion and the
output indicated it failed.  Very often tests are reported to fail because of a
`Timeout`, which means the test ran, but longer than was allowed.  If
this happens, the test timeout needs to be increased. When configuring GridPACK,
add an option to cmake like 

```
    -D GRIDPACK_TEST_TIMEOUT:STRING=30
```

which indicates any unit is only allowed to run for 30 seconds before is killed.
Another common failure is reported by `ctest` as `Required
regular expression not found`. This means that the test ran to completion,
but the output was not recognized as a correct.  Most of the time this indicates
a real error, but sometimes, particularly with parallel tests, 

If lots of tests fail, there is probably some problem with the GridPACK
configuration or some other system setup (like MPI peculiarities). If only a few
tests fail, they can probably be ignored, unless the test failure indicates some
functionality of importance to you.  If a test fails, it possible that a bug
report has been made. Check [current
issues](https://github.com/GridOPTICS/GridPACK/issues) to see. 

Failed tests can be diagnosed by running them individually and looking at the
output.  For example, to run test 2 only, run `ctest` at the top of
the GridPACK build directory like this

```
    ctest -VV -I 2,2
```

This produces

```
    r/src/build/DartConfiguration.tcl
    UpdateCTestConfiguration  from
:/home/d3g096/Projects/GridPakLDRD/GridPACK.github/master/src/build/DartConfiguration.tcl
    Test project
/home/d3g096/Projects/GridPakLDRD/GridPACK.github/master/src/build
    Constructing a list of tests
    Done constructing a list of tests
    Checking test dependency graph...
    Checking test dependency graph end
    test 2
        Start 2: greetings_parallel
    
    2: Test command: /usr/lib64/openmpi/bin/mpiexec "-n" "4" "greetings"
    2: Test timeout computed to be: 9.99988e+06
    2: I am process 0 of 4.
    2: 
    2: Creating communicators using split:
    2: 
    2: I am process 1 of 4.
    2: I am process 2 of 4.
    2: I am process 3 of 4.
    2: I am process 0 (original process is 0) of 2 on communicator 0.
    2: I am process 1 (original process is 1) of 2 on communicator 0.
    2: I am process 0 (original process is 2) of 2 on communicator 1.
    2: I am process 1 (original process is 3) of 2 on communicator 1.
    2: 
    2: Creating communicators using divide:
    2: 
    2: I am process 0 (original process is 2) of 2.
    2: (assignment) I am process 0 (original process is 2) of 2.
    2: (copy) I am process 0 (original process is 2) of 2.
    2: I am process 1 (original process is 3) of 2.
    2: (assignment) I am process 1 (original process is 3) of 2.
    2: (copy) I am process 1 (original process is 3) of 2.
    2: I am process 1 (original process is 1) of 2.
    2: (assignment) I am process 1 (original process is 1) of 2.
    2: (copy) I am process 1 (original process is 1) of 2.
    2: I am process 0 (original process is 0) of 2.
    2: 
    2: Testing assignment operator for communicators:
    2: 
    2: (assignment) I am process 0 (original process is 0) of 2.
    2: 
    2: Testing copy constructor for communicators:
    2: 
    2: (copy) I am process 0 (original process is 0) of 2.
    2: 
    2: Testing summation operator for communicators
    2: 
    2: 
    2: Summation test for communicators passed
    2: 
    1/1 Test #2: greetings_parallel ...............   Passed    0.31 sec
    
    100% tests passed, 0 tests failed out of 1
    
    Total Test time (real) =   0.32 sec
```

If you wish to report a bug about a failed unit test, run the test as described
and post it with the bug report.  

## Running the Powerflow Example(s)

If desired, the powerflow example application can be run on a number of
different networks.  Starting from the top GridPACK build directory, change into
the powerflow application build directory:

```
    cd applications/powerflow
```

Run the IEEE 14 bus problem serially using

```
    ./pf.x input_14.xml
```

or in parallel using

```
    mpiexec -np 2 ./pf.x input_14.xml
```

(but don't use too many processors).  The important part of the output from this
case is

```
   Branch Power Flow
   
        Bus 1       Bus 2   CKT         P                    Q
          1           2     BL     156.882891           -20.404292
          1           5     BL      75.510382             3.854991
          2           3     BL      73.237579             3.560203
          2           4     BL      56.131496            -1.550350
          2           5     BL      41.516215             1.170998
          3           4     BL     -23.285690             4.473116
          4           5     BL     -61.158230            15.823642
          4           7     BL      28.074176            -9.681066
          4           9     BL      16.079758            -0.427611
          5           6     BL      44.087321            12.470680
          6          11     BL       7.353277             3.560473
          6          12     BL       7.786067             2.503414
          6          13     BL      17.747977             7.216575
          7           8     BL      -0.000000           -17.162971
          7           9     BL      28.074176             5.778691
          9          10     BL       5.227552             4.219138
          9          14     BL       9.426381             3.610006
         10          11     BL      -3.785322            -1.615063
         12          13     BL       1.614258             0.753959
         13          14     BL       5.643851             1.747174
   
   Generator Power
   
   Bus Number  GenID        Pgen              Qgen
          1       1       2.323933         -0.165493
          2       1       0.400000          0.435571
          3       1       0.000000          0.250753
          6       1       0.000000          0.127309
          8       1       0.000000          0.176235
   
   Bus Voltages and Phase Angles
   
   Bus Number      Phase Angle      Voltage Magnitude
          1          0.000000             1.060000
          2         -4.982589             1.045000
          3        -12.725100             1.010000
          4        -10.312901             1.017671
          5         -8.773854             1.019514
          6        -14.220946             1.070000
          7        -13.359627             1.061520
          8        -13.359627             1.090000
          9        -14.938521             1.055932
         10        -15.097288             1.050985
         11        -14.790622             1.056907
         12        -15.075585             1.055189
         13        -15.156276             1.050382
         14        -16.033645             1.035530
```

Input for other networks are available:

* `input_118.xml`: a (IEEE?) 118 bus case
* `input_polish.xml`: Polish network model with 3120 buses and 3684 branches
* `input_european.xml`: European network model with 13659 buses and 18625 branches

for which more processors can be deployed.