## PETSc

GridPACK currently relies on the [Portable, Extensible Toolkit for Scientific
Computation (PETSc)][http://www.mcs.anl.gov/petsc/index.html) for parallel
linear algebra, and linear and nonlinear system solvers. The PETSc interface
tends to change a bit as new releases come out, requiring adjustments in any
applications that use it. We currently recommend using PETSc version 3.16 with
GridPACK.

PETSc is a complicated package with numerous options.  PETSc needs to be built
with MPI enabled and using the same MPI implementation used for GridPACK and the
other libraries such as Boost and GA.
Originally, GridPACK could only use PETSc if it was configured for
complex support. The current GridPACK release can use either complex or real
builds. However, most applications in GridPACK use complex matrices, so it is
still preferable to configure PETSc to use complex variables.
Refer to the
[PETSc installation documentation](http://www.mcs.anl.gov/petsc/documentation/installation.html)
for additional information on how to configure PETSc.

Configuring and building PETSc is done in the top level PETSc directory. One of
the configuration variables that needs to be set when configuring and building
PETSc is PETSC_ARCH. In the example below, PETSC_ARCH was set to
`'arch-linux2-complex-opt'`. After the build is complete, there will be a
directory beneath the top level directory with whatever name was assigned to
PETSC_ARCH. This directory contains the include and lib directories for the
PETSc libraries.

The GridPACK configuration must know where
[PETSc](http://www.mcs.anl.gov/petsc/index.html) is installed.  This is specified
by two options as shown below. 

    -D PETSC_DIR:STRING='/Users/d3g096/ProjectStuff/petsc-3.16.3' \
    -D PETSC_ARCH:STRING='arch-linux2-complex-opt' \

Currently, the configuration will recognize and adjust the GridPACK build if the
[PETSc](http://www.mcs.anl.gov/petsc/index.html) build includes
[ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview),
[SuperLU_DIST](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) and/or
[MUMPS](http://mumps.enseeiht.fr/).  Many of the example GridPACK applications
expect a parallel direct linear solver to be built into PETSc.  This is
satisfied by including
[SuperLU_DIST](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/)
or [MUMPS](http://mumps.enseeiht.fr/) in the PETSc build.

A typical configuration script for building PETSc is

```
    ./configure \
       PETSC_ARCH=arch-linux2-complex-opt \
       --with-scalar-type=complex \
       --download-superlu_dist \
       --download-superlu \
       --download-parmetis \
       --download-metis \
       --download-suitesparse \
       --download-f2cblaslapack \
       --with-clanguage=c++ \
       --with-shared-libraries=0 \
       --with-x=0 \
       --with-mpiexec=mpiexec \
       --with-debugging=0
```

The PETSC_ARCH variable can be set to whatever the user wants but should be
descriptive of this particular build. PETSc can be reconfigured with another
value of the PETSC_ARCH variable to support multiple different PETSc builds.

At the end of the the configuration, you should see something similar to the
following output

```
xxx=========================================================================xxx
 Configure stage complete. Now build PETSc libraries with:
   make PETSC_DIR=/pic/projects/gridpack/software/petsc-3.16.3 PETSC_ARCH=linux-openmpi-gnu-cxx-complex-opt all
xxx=========================================================================xxx
```

Cut and paste the `make` command into the command line to compile the PETSc
libraries. At the end of the build you will see

```
=========================================
Now to check if the libraries are working do:
make PETSC_DIR=/pic/projects/gridpack/software/petsc-3.16.3 PETSC_ARCH=linux-openmpi-gnu-cxx-complex-opt check
=========================================
```

At this point, the build is done, but if you want to test the libraries, you can
cut and paste the make command to run some quick tests to verify
that the build was successful.

If you want to build PETSc with shared libraries, change the argument of
`--with-shared-libraries=` from 0 to 1. To build PETSc using
real variables instead of complex, set the argument of `--with-scalar-type` from
`complex` to `real`.

## Cross-Compilation

On some platforms, the backend hardware is different from the front end, meaning
that software compiled on the front end can only be run on the back end (compute
nodes). This is a problem for the PETSc build, which creates small programs and
runs them as part of the configuration procedure to see identify properties of
the system.

The work-around for this is to configure PETSc using the

```
--with-batch=1
```
option. In older versions of PETSc, this would cause the configure process to
generate a script with a name based on the `PETSC_ARCH` variable, (e.g.
`conftest-arch-linux2-complex-opt`), that would then be submitted to the batch
queuing system. After submission, the configure would then proceed normally. For
newer versions of PETSc, including 3.16.x, this does not appear to be necessary,
and including the `--with-batch` flag is the only modification necessary when
cross-compiling.
