Installation
============

The CLI is the ``gridpack`` console script installed with the ``gridpack``
Python package, so installing the package installs the CLI.  These are the
package's steps; ``python/README.md`` in the source tree is the reference.

Prerequisites
-------------

- GridPACK built and *installed* as shared libraries, with PETSc, Global
  Arrays and Boost also shared
- Python 3.9 to 3.13, ``pip``, ``setuptools`` >= 61
- ``mpi4py`` built against the same MPI as GridPACK
- The compiler, CMake and MPI used to build GridPACK
- The ``pybind11`` submodule: ``git submodule update --init``

Install
-------

Point the build at the GridPACK installation and install into a virtual
environment::

    export GRIDPACK_DIR=/path/to/gridpack/install   # holds lib/GridPACK.cmake
    python3 -m venv ~/.venvs/pygridpack
    source ~/.venvs/pygridpack/bin/activate
    pip install setuptools wheel mpi4py
    cd /path/to/GridPACK
    pip install --no-build-isolation -e python/

``GRIDPACK_DIR`` is not optional: unset, the build falls back to a
hard-coded developer path and fails confusingly.  ``--no-build-isolation``
makes the build use the ``mpi4py`` already in the environment rather than
compiling a fresh one against whatever MPI ``pip`` finds first.

Verify
------

::

    gridpack --version
    gridpack powerflow input_14.xml
    cd python && pytest -q

Every test should pass and none should skip: a skip means ``mpiexec`` is
not on ``PATH`` or the extension did not build.  ``pytest --fail-on-skip``
turns a skip into a failure; the Docker build uses it as its gate.

For the published image with the whole stack already built, see
:doc:`docker`.
