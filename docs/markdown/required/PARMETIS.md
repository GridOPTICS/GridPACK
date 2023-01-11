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

In the event that you cannot build ParMETIS as part of PETSc, you can configure
and build ParMETIS separately (this should rarely, if ever, be necessary).

It's easiest to include ParMETIS in the PETSc build. The GridPACK configuration
will recognize and use ParMETIS from the PETSc installation. However, if you
want to build ParMETIS separately, instructions are below. Information on
downloading tar files and creating scripts can be found
[here](LINUX_BASICS.md$linux-basics).

Starting in the ParMETIS source directory, build and install METIS first:

```
    cd metis
    make config prefix="$PREFIX"
    make
    make install
```

The variable `$PREFIX` should be replaced with location where you want to
install the METIS and PARMETIS libraries. Then build and install ParMETIS:

```
    cd ..
    make config cc=mpicc cxx=mpicxx prefix="$PREFIX"
    make 
    make install
```

Do some tests to make sure it works:

```
    cd Graphs
    mpirun -np 4 ptest rotor.graph rotor.graph.xyz
    mpirun -np 4 ptest rotor.graph
    mpirun -np 8 ptest bricks.hex3d
```

The last test appears to hang but does not to otherwise cause any
problems.

In order to get ParMETIS 4.0 to compile with older GNU compilers, a warning
option needs to be removed from one of the build system files.  In the top
ParMETIS source directory, execute the following command:

```
    sed -i.org -e 's/-Wno-unused-but-set-variable//g' metis/GKlib/GKlibSystem.cmake
```
