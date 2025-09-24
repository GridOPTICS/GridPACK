Configuring and Building GridPACK
=================================

**A note about CMake**: The command for invoking CMake in this manual
and the documentation in https://github.com/GridOPTICS/GridPACK is
usually of the form

::

   cmake [OPTIONS] ..

This particular form assumes that the build directory is below the
directory that contains the top-level ``CMakeLists.txt`` file for
the build. For GridPACK, this is located in the ``src`` directory.
If your build directory for GridPACK is below ``src`` and you invoke
CMake from this directory, the “``..``” at the end of the
``cmake`` command is pointing to ``src``. You could also use the
absolute path to the ``src`` directory instead of “``..``” and
this would work no matter where you locate the build directory.

Build Scripts
-------------

Most of the documentation for building GridPACK describes how to build
individual libraries used by GridPACK and how to link these to the GridPACK
executables. See the section on Manual Builds below. Recently, we have
provided scripts that automate much of this process. These scripts work on many
platforms but are not guaranteed to work all the time. In the event of a failure
you will either have to see if you can get the scripts working by fixing
whatever bugs you may encounter or by reverting to doing a manual build as
decribed below.

The two scripts used to build GridPACK are located in the top level directory
and are called ``install_gridpack_deps.sh`` and ``install_gridpack_deps.sh``.
These scripts can be run by simply typing

::

   install_gridpack_deps.sh

and

::

   install_gridpack.sh

in the top level directory.
The first script will download and build all the libraries that GridPACK depends
on. These include Boost, PETSc and Global Arrays (GA). This does not include
MPI. Most HPC clusters will have MPI already installed but for laptops and
workstations, you may need to build or install these libraries yourself. If the
``install_gridpack_deps.sh`` script is used to install the external
dependencies, the libraries will be located in a directory call
``external-dependencies`` that is located under the top level GridPACK
directory. Note that this script may take around an hour to run, depending on
the system.

Once the external dependences have been built, the
GridPACK framework and its associated applications are built by running the
``install_gridpack.sh`` script. If you want to include the python interface,
then will need to edit the top of this script and set both the
``install_gridpack_shared`` and ``install_gridpack_python1`` variables to true.
Once this script has been run, two directories shou appear under the
``GRIDPACK/src`` directory. The first is ``build``, which will contain all the
application executables in GridPACK and the second is ``install``. The install
directory is chiefly of interest for users that are interested in developing
their own applications and contains the libraries and header files for
individual framework components. These can be used to compile and link new
GridPACK-based applications.

These scripts build the Global Arrays library with the two-sided runtime. This
is straightforward to use but is low performing on large numbers of processors.
Better performance can be achieved with the progress ranks runtime. To build
with this runtime, edit the ``install_gridpack_deps.sh`` file and replace the
string ``--with-mpi-ts`` with ``-with-mpi-pr`` in the section on building GA. In
addition, you will need to add the line

::

   -D USE_PROGRESS_RANKS:BOOL=TRUE

to the list of ``cmake_args`` in ``install_gridpack.sh``. For additional
information on running applications with this runtime, see the information on
using progress ranks in the section below.

Manual Builds
-------------

Building GridPACK requires several external libraries that must be built
prior trying to configure and build GridPACK itself. On some systems,
these libraries may already be available but in many cases, users will
need to build them by hand. An exception is MPI, which is usually
available on parallel platforms, although users interested in running
parallel jobs on a multi-core workstation may still need to build it
themselves. In any case, the best way to guarantee that all libraries
are compatible with each other is to build them all using a consistent
environment and set of compilers. There is extensive documentation on
how to build GridPACK and the libraries on which it depends on the
website located at https://github.com/GridOPTICS/GridPACK. We refer to
the information on the website for most of the details on how to build
GridPACK and will only discuss some general properties of the configure
procedure in this document.

Example scripts for building the libraries used by GridPACK on different
systems can be found under ``$GRIDPACK/src/scripts``. In most cases
these need to be modified slightly before they will work on your system,
but the changes are usually small and self-evident. The scripts contain
some additional documentation at the top to help you with these
modifications. Find a script for a platform that is similar to your
system and use this as the starting point for your build. Addititional
information for building on advanced platforms, such as DOE’s Leadership
Class Facilities, can be found in the ``$GRIDPACK/docs/notes``
directory. This directory contains notes on building GridPACK on
different platforms that may be useful to users trying to build on
similar systems. Note that these systems can be quite complicated to
build on for any application and usually require assistance from
facility staff.

GridPACK uses the CMake build system to create a set of make files that
can then be used to compile the entire GridPACK framework. Most of the
effort in building GridPACK is focused on getting the configure process
to work, once configure has been successfully completed, compilation is
usually straightforward. Builds of GridPACK should be done in their own
directory and this also makes it possible to have multiple builds that
use different configuration parameters associated with the same source
tree. Typically, the build directories are under ``$GRIDPACK/src``
directory but they can be put anywhere the user chooses. The user then
needs to run CMake from the build directory to configure GridPACK and
then ``make`` and ``make install`` to compile and install the
GridPACK libraries. After running ``make``, all applications in the
GridPACK source tree are also available for use. The application
executables will be located in the build directory and not in the source
tree.

GridPACK currently makes use of five different libraries. MPI and Global
Arrays are used for communication, Boost provides several C++ extensions
used throughout GridPACK, Parmetis is used to partition networks over
multiple processors and PETSc provides parallel solvers and algebraic
functionality. Except for MPI, which is usually available through
compiler wrappers such as ``mpicc`` and ``mpicxx``, the
locations of the remaining libraries need to be specified in the CMake
configure command.

Because the cmake command takes a large number of arguments, it is
usually a good idea to put the entire command in a script. The script
can then be edited as needed. Make sure that the script is executable by
running the ``chmod +x`` command on it. A typical CMake configure
script is

::

   rm -rf CMake*

   cmake -Wdev \
         -D BOOST_ROOT:STRING='$HOME/software/boost_1_65_0' \
         -D PETSC_DIR:STRING='$HOME/software/petsc-3.16.0' \
         -D PETSC_ARCH:STRING='linux-openmpi-gnu-cxx' \
         -D GA_DIR:STRING='$HOME/software/ga-5-8' \
         -D USE_PROGRESS_RANKS:BOOL=FALSE \
         -D MPI_CXX_COMPILER:STRING='mpicxx' \
         -D MPI_C_COMPILER:STRING='mpicc' \
         -D MPIEXEC:STRING='mpiexec' \
         -D CMAKE_INSTALL_PREFIX:PATH='$GRIDPACK/src/build/install' \
         -D CMAKE_BUILD_TYPE:STRING='RELWITHDEBINFO' \
         -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
         -D CMAKE_VERBOSE_MAKEFILE:STRING=TRUE \
         ..

The first line removes any configuration files that may be left over
from a previous configuration attempt. Removing these files is generally
a good idea since parameters from a previous unsuccessful attempt may
bleed over into the current configuration and either spoil the
configuration itself or lead to problems when you try to compile the
code.

The Boost, PETSc, Parmetis and Global Array library locations are
specified by the ``BOOST_ROOT``, ``PETSC_DIR`` and
``GA_DIR`` variables. In most cases, Parmetis can be built with
PETSc and will be detected by the CMake build system. In the rare case
that this is not true, the location of Parmetis can be specified with
the ``PARMETIS_DIR`` variable. The ``PETSC_ARCH`` variable
specifies the particular build within PETSc that you want GridPACK to
use.

The Global Arrays library can be built using a number of different
runtimes. The default runtime uses MPI two-sided communication. While it
is very easy to use, this runtime does not scale well beyond a dozen or
so processors. Users interested on running on large numbers of cores
should look at configuring Global Arrays with other runtimes. A high
performing GA runtime that is available on all platforms is called
progress ranks. While it is very robust, this runtime has a peculiarity
in that it reserves one MPI process per SMP node to manage
communication. Thus, if you request a total of 20 MPI processes on 4
nodes with 5 processes running on each node only 4 MPI process per node
will actually be available to the application for a total of 16. In
order to notify GridPACK that you are using this runtime, you need to
set the parameter ``USE_PROGRESS_RANKS`` to true. In the example
above, we are not using progress ranks so we set
``USE_PROGRESS_RANKS`` to false. The ``GA_EXTRA_LIBS`` parameter
can be used to include extra libraries in the link step that are not
picked up as part of the configuration process. For this example, this
parameter is not needed.

The MPI wrappers for the C and C++ compilers can be specified by setting
``MPI_C_COMPILER`` and ``MPI_CXX_COMPILER`` and the MPI launch
command can be specified using ``MPIEXEC``. The
``CMAKE_INSTALL_PREFIX`` specifies the location of the installed
build of GridPACK. This location should be used when linking external
applications to GridPACK. The ``CMAKE_BUILD_TYPE`` can be used to
control the level of debugging symbols in the library.
``MPIEXEC_NUM_PROCS`` should be set to a small number and controls
the number of processors that will be used if running the parallel tests
in the GridPACK test suite. Many of the application tests are small (9
or 14 buses) and will fail if you try and run on a large number of
cores. Finally, ``CMAKE_VERBOSE_MAKEFILE`` controls the level of
information generated during the compilation. It is mainly of interest
for people doing development in GridPACK and most other users can safely
set it to false.

A new feature in the build is to use shared libraries instead of static
builds. This may be of interest to users that are interested in wrapping
GridPACK applications with python. A shared library version of GridPACK
can be created by configuring GridPACK against versions of Boost, GA,
and PETSc that are built as shared libraries. It appears that just
configuring against shared libraries is enough to trigger a share
library build in CMake, but users can add the line

::

   -D BUILD_SHARED_LIBS:BOOL=TRUE \

to their configuration invocation to make sure.

The final argument of the cmake command is the location of the top level
``CMakeLists.txt`` file in the source tree. For GridPACK, this file
is located in the ``$GRIDPACK/src`` directory. The above example
assumes that the build directory is located directly under
``$GRIDPACK/src`` so the ``..`` at the end of the configure
script is pointing to the directory containing the
``CMakeLists.txt`` file.

Running GridPACK Applications
-----------------------------

Once the GridPACK framework has been built, applications and framework
tests can be run using standard MPI scripts for running jobs. A typical
invocation to run a code ``code.x`` on some number of processors is

::

   mpirun -n 2 code.x input.xml

The command ``mpirun`` is used to launch an application on multiple
processors and invokes the MPI library. The ``-n`` indicates the number
of processors that the run will use.  In this case the code will run on
2 processors.  Note that different MPI implementations may use different
commands for launching MPI jobs.  Another common command is ``mpiexec``.
Consult your local system documentation for details. Applications may
also have additional arguments that are processed inside the application
itself. Most GridPACK applications will take an argument representing
the input file for the application. In this example, the input file is
``input.xml``. Additional information needed to execute the job is
typically located within the input XML file, so GridPACK application do
not take more than one argument.
