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

* [Redhat Linux Workstation](../DUMMY.md)
