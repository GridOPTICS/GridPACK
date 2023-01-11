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

* [Mac Yosemite ](../DUMMY.md)
* [Mac High Sierra](../DUMMY.md)
* [Redhat Linux Workstation](../DUMMY.md)
* [CentOS 6](../DUMMY.md)
* [Redhat Linux Cluster](../RC_CLUSTER.md)
