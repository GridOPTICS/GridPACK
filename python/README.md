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

It is *not* necessary with MPICH or apparently OpenMPI on other Linux
distributions.  

### Build

*This is deprecated, but may still work with older Pythons.*

In the top directory (where `setup.py` is),
```
python setup.py build
```

### Test

*This is deprecated, but may still work with older Pythons.*

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

With modern versions of Python and setuptools, it's best to use `pip`
to install the module:
```
setenv GRIDPACK_DIR $prefix
cd .../GridPACK/python
pip install --no-deps --upgrade --prefix=$GRIDPACK_DIR .
```
This will install the Python module in the same place as the rest of
GridPACK (`$GRIDPACK_DIR`).  To use the module, `PYTHONPATH` must be
set to the site packages directory under `$GRIDPACK_DIR`. The path
will depend on the Python version used.  For example, 
```
setenv PYTHONPATH $GRIDPACK_DIR/lib/python3.12/site-packages/
```
or
```
export PYTHONPATH=$GRIDPACK_DIR/lib/python3.12/site-packages
```
Alternatively, a Python virtual environment can be used.  

An easy check to see if the module works is
```
python -c 'import gridpack'
```

*The below deprecated, but may still work with older Python
installations. If this is used with newer Python versions, the module
will not be installed correctly and be unusable.  *

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

## Run Examples

Once installed, with `PYTHONPATH` set correctly, the Python version of
several GridPACK applications can be run.  In
`.../python/src/example` are a number of custom Python scripts that
can be run to insure the viability of the GridPACK Python
interface:
```
cd src/example
python 39bus_test_example.py
python 39bus_test_example_dsf.py
python 39bus_scatterload_steptest_new_itr.py
mpiexec -np 2 python 39bus_scatterload_steptest_new_itr.py
python 39bus_scatterload_steptest_new_itr_dsf.py
mpiexec -np 2 python 39bus_scatterload_steptest_new_itr_dsf.py
python 39bus_scatterload_steptest_new_itr_compensateY.py
python 39bus_test_pfdata.py
```
Dynamic simulation examples can also be run.
`$GRIDPACK_DIR/bin/dsf2.py` behaves identically to
`$GRIDPACK_DIR/bin/dsf2.x` and can run any of the example simulations
in `.../src/build/applications/dynamic_simulation_full_y`, for example

```
cd ../src/build/applications/dynamic_simulation_full_y
$GRIDPACK_DIR/bin/dsf2.py input_9b3g.xml
$GRIDPACK_DIR/bin/dsf2.py input_240bus.xml
mpiexec -np 4 $GRIDPACK_DIR/bin/dsf2.py input_240bus.xml
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


