Docker
======

``Dockerfile.pygridpack`` at the repository root adds the Python package
and the CLI to the published ``pnnl/gridpack`` image, which already holds
GridPACK, PETSc, Global Arrays, Boost and MPI.  It builds in minutes rather
than the hours the full ``Dockerfile`` takes::

    docker build -f Dockerfile.pygridpack -t pygridpack .
    docker run --rm -it -v "$PWD:/app/workspace" pygridpack

Inside the container no environment variables are needed::

    gridpack powerflow input_14.xml
    mpiexec -np 4 gridpack powerflow input_14.xml

The build refuses to produce an image whose import, ``gridpack --version``
or power flow smoke fails, then runs the full test suite with
``--fail-on-skip``, so a skipped test fails the build as well.  The suite
can also be run against a finished image::

    docker run --rm pygridpack pytest /app/python/tests -q

Inputs and outputs live in the mounted directory: XML, RAW and DYR files
go in, and what GridPACK writes (power flow tables, contingency ``.out``
files, generator watch CSVs) comes out there.  ``mpiexec`` inside the
container sets the rank count; Docker replicas do not.
