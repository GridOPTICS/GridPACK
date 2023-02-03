The following packages are required in order to build GridPACK. Some of them may
already be available on your platforms, either as modules or as part of the
system software. Be aware that in some cases, the software may be too old or may
not have been built with features required by GridPACK. Existing CMake
executable may be too old and installations of Boost are frequently missing the
Boost MPI libraries. If this is the case, then it may be necessary to upgrade
the existing software.

The following functionality is require by GridPACK

- [CMake](LINUX_BASICS.md#cmakectest) Gridpack currently requires version 3.5.0
  or greater.
- [MPI](LINUX_BASICS.md#mpi) Any recent version of MPI should work. Most of the
  builds described in these pages were done using
  [OpenMPI](http://www.open-mpi.org/)
- [Boost](BOOST.md) Note that you will need a Boost library built with Boost
  MPI. Many preinstalled versions of Boost and Boost modules do not contain this
  library. If they don't, you will need to build your own version of Boost.
- [Global Arrays](GLOBAL_ARRAYS.md) Any recent version of Global Arrays should
  work. The current release is 5.8.
- [PETSc](PETSC.md) There have been some substantial changes in how PETSc
  supports builds since GridPACK was originally written. We are currently
  supporting more recent versions of PETSc (3.16) and are deprecating the use of
  older versions.
- [ParMETIS](PARMETIS.md) This library can, in most circumstances, be included
  as part of the PETSc libraries. In the rare cases that this is not possible,
  ParMETIS needs to be built and included separately.
