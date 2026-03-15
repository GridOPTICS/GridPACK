Installation
============

This chapter describes how to build, install, and verify the GridPACK
CLI and its Python bindings.

Prerequisites
-------------

The following software must be installed before building:

- **Python** 3.6 or later
- **MPI** implementation (MPICH, OpenMPI, or equivalent)
- **mpi4py** Python package
- **GridPACK** source code with external dependencies:

  - PETSc (with SuperLU_DIST)
  - GA (Global Arrays) 5.9 or later, built with shared libraries
  - Boost (MPI, serialization)

Building GridPACK
-----------------

GridPACK should be built with shared libraries enabled to ensure
compatibility between the C++ executables and the Python bindings.

Step 1: Build GA with Shared Libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If GA was built with ``--disable-shared``, rebuild it::

    cd /path/to/external-dependencies/ga-5.9
    make clean
    ./configure --prefix=$(pwd)/install_for_gridpack \
        --with-mpi-ts --disable-f77 --without-blas \
        --enable-cxx --enable-i4 --enable-autodetect \
        --with-pic --enable-shared --enable-static \
        CC=mpicc CXX=mpicxx
    make -j4
    make install

Step 2: Configure and Build GridPACK
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configure GridPACK with shared libraries and GA 5.9::

    cd /path/to/GridPACK/src/build
    cmake .. \
        -DBUILD_SHARED_LIBS=ON \
        -DGA_DIR=/path/to/external-dependencies/ga-5.9/install_for_gridpack \
        -DGA_LIBRARY=/path/to/external-dependencies/ga-5.9/install_for_gridpack/lib/libga.dylib \
        -DGA_CXX_LIBRARY=/path/to/external-dependencies/ga-5.9/install_for_gridpack/lib/libga++.dylib \
        -DARMCI_LIBRARY=/path/to/external-dependencies/ga-5.9/install_for_gridpack/lib/libarmci.dylib \
        -DGA_INCLUDE_DIR=/path/to/external-dependencies/ga-5.9/install_for_gridpack/include
    make -j4
    make install

Step 3: Build the Python Bindings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Python extension module is built separately from the main
GridPACK build::

    cd /path/to/GridPACK/python/build_manual
    cmake .. \
        -DPython_EXECUTABLE=$(which python3) \
        -DGRIDPACK_DIR=/path/to/GridPACK/src/install \
        -DMPI_CXX_INCLUDE_DIRS=/path/to/mpi/include
    make -j4

Then install the compiled module::

    cp src/gridpack.cpython-*.so \
        /path/to/GridPACK/src/install/lib/python3.X/site-packages/

macOS Code Signing
~~~~~~~~~~~~~~~~~~~

On macOS, newly built shared libraries must be ad-hoc code signed
or the operating system will terminate any process that tries to
load them. Sign all GridPACK and GA shared libraries after building::

    # Sign GA libraries
    codesign --force --sign - /path/to/ga-5.9/install_for_gridpack/lib/libga.1.dylib
    codesign --force --sign - /path/to/ga-5.9/install_for_gridpack/lib/libga++.0.dylib
    codesign --force --sign - /path/to/ga-5.9/install_for_gridpack/lib/libarmci.0.dylib
    codesign --force --sign - /path/to/ga-5.9/install_for_gridpack/lib/libcomex.0.dylib

    # Sign GridPACK libraries
    for f in /path/to/GridPACK/src/install/lib/libgridpack_*.dylib; do
        codesign --force --sign - "$f"
    done

    # Sign the Python extension module
    codesign --force --sign - /path/to/site-packages/gridpack.cpython-*.so

.. note::

   If ``import gridpack`` is silently killed (exit code 137 / signal 9)
   without any error message, unsigned shared libraries are almost
   certainly the cause. Re-run the signing commands above.

Step 4: Verify the Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set the Python path and test the import::

    export PYTHONPATH=/path/to/GridPACK/src/install/lib/python3.X/site-packages
    python3 -c "
    import gridpack
    print('Submodules:', [m for m in dir(gridpack) if not m.startswith('_') and m[0].islower()])
    c = gridpack.powerflow.Contingency()
    print('Contingency class: OK')
    "

Expected output::

    Submodules: ['dynamic_simulation', 'emt', 'hadrec', 'powerflow', 'state_estimation']
    Contingency class: OK

To test the CLI::

    export PYTHONPATH=/path/to/site-packages:/path/to/GridPACK/python/src
    python3 /path/to/GridPACK/python/src/gridpack_cli.py --help

Installing as a Console Script
-------------------------------

To make the ``gridpack`` command available system-wide without
setting ``PYTHONPATH`` manually, install the Python package::

    cd /path/to/GridPACK/python
    pip install .

This registers the ``gridpack`` console script entry point. After
installation, the command is available directly::

    gridpack --help
    gridpack powerflow input.xml

For development (editable) installation::

    pip install -e .

Rebuilding After Source Changes
-------------------------------

If you modify the pybind11 binding source files (e.g.,
``python/src/gridpack_pf.cpp``), rebuild and reinstall the Python
extension module::

    cd /path/to/GridPACK/python/build_manual
    make -j4
    cp src/gridpack.cpython-*.so /path/to/site-packages/
    codesign --force --sign - /path/to/site-packages/gridpack.cpython-*.so

Environment Configuration
--------------------------

If GridPACK is installed to a non-standard location, ensure that
the shared libraries and Python module are discoverable::

    export PYTHONPATH=/path/to/GridPACK/src/install/lib/python3.X/site-packages:/path/to/GridPACK/python/src

On Linux::

    export LD_LIBRARY_PATH=/path/to/GridPACK/src/install/lib:$LD_LIBRARY_PATH

On macOS, ``DYLD_LIBRARY_PATH`` may be used but is generally not
needed if the libraries are built with correct ``rpath`` settings
(which is the default for the CMake build).
