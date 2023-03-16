# GridPACK Python Wrapper

This provides Python wrappers to a limited set of the
[GridPACK](https://www.gridpack.org) library.  The HADREC application
is the main capability currently exposed.

## Requirements

  * GridPACK >= 3.4
  * Python >= 2.7 (but use >= 3.x)
  * pybind11 >= 2.4
  * Python `setuptools` package
  * Python `nose` package
  * CMake used to build GridPACK
  * C++ Compiler used to build GridPACK
  * MPI (and other library) installation used to build GridPACK

### Pybind11

The source for `pybind11` is needed to build the wrapper.
Consequently, it's source is included as a submodule.  To make it's
there before building update the submodule:
```
git submodule update --init
```

### GridPACK

GridPACK must be built and *installed* as *shared* libraries. This
requires that any GridPACK dependencies (e.g. PETSc, Global Arrays,
Boost) also be built as shared libraries.

## Build and Test

### Environment Variables

There are a couple of environment variables that control the Python
wrapper build. 

First, the `GRIDPACK_DIR` environment variable is used to indicate
where GridPACK was installed.  So, in Bourne-like shells,
```
GRIDPACK_DIR=/usr/local/gridpack
export GRIDPACK_DIR
```

Second, on RHEL systems (probably CentOS too), the stock OpenMPI
installation causes problems.  To alleviate these set
`RHEL_OPENMPI_HACK` to anything, e.g.

```
RHEL_OPENMPI_HACK=yes
export RHEL_OPENMPI_HACK
```

This seems to be necessary on constance.

### Build

In the top directory (where `setup.py` is),
```
python setup.py build
```

### Test
Again, in the top directory,
```
python setup.py test
```
This will run 3 tests.  One test is the full HADREC application on the
TAMU network.  This test spews lots of text to the terminal, but
should end with something like
```
ok
hello_test (tests.gridpack_test.GridPACKTester) ... ok
task_test (tests.gridpack_test.GridPACKTester) ... ok

----------------------------------------------------------------------
Ran 3 tests in 7.547s

OK
```

## Installation 

In order to use the wrapper it needs to be installed in a way Python
can use it.  However, you shouldn't (and probably cannot) install the
module in the system python library.  So, choose a place to install it. I would
recommend the same directory as GridPACK (i.e. `$GRIDPACK_DIR` as used
above).  

The installation process forces the `PYTHONPATH` environment to be
correct *before* installation. If `PYTHONPATH` is not empty 
```
PYTHONPATH="${GRIDPACK_DIR}/lib/python:${PYTHONPATH}"
export PYTHONPATH
python setup.py install --home="$GRIDPACK_DIR"
```
This seems to be the way on Linux systems. 

This installs the Python module and a Python version of the HADREC
application, `${GRIDPACK_DIR}/bin/hadrec.py`.

## Run Example

Once installed, with `PYTHONPATH` set correctly, the Python version of
the HADREC application can be run

```
cd src/tests
python ${GRIDPACK_DIR}/bin/hadrec.py input_tamu500_step005.xml
```

## Running on Constance

I've built GridPACK and the HADREC wrapper on constance.  

Use the following modules:
```
module purge
module load gcc/4.9.2
module load openmpi/1.8.3
module load python/2.7.3
module load cmake/3.8.2
```
Use the following environment variables:
```
GRIDPACK_DIR=/pic/projects/gripdack/hadrec/gridpack
PYTHONPATH="$GRIDPACK_DIR/lib/python"
PATH="${GRIDPACK_DIR}/bin:$PATH"
```
With these settings, the Python version of
the HADREC application can be run with
```
cd src/tests
hadrec.py input_tamu500_step005.xml
```


