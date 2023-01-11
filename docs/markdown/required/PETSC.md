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

* [Mac Yosemite ](../DUMMY.md)
* [Mac High Sierra](../DUMMY.md)
* [Redhat Linux Workstation](../DUMMY.md)
* [CentOS 6](../DUMMY.md)
* [Redhat Linux Cluster](../RC_CLUSTER.md)
