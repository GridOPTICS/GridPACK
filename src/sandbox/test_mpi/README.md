# Simple PyBind11 Wrapper

## Requirements

  * Python >= 2.7
  * pybind11 >= 2.4
  * Python `setuptools` package
  * Python `nose` package
  * CMake
  * C++ Compiler
  * MPI 

If pybind11 is installed, but without CMake support, the build will fail.  In that event, place the pybind11 code in the top directory:
```
git clone -b v2.4 https://github.com/pybind/pybind11.git
```
It is not needed after the module is built.
  
## Build and Test

Build the module:
```
python setup.py build
```
Test the module by running test (need `nose` for this):

```
python setup.py test
```
This should produce output like
```
test_is_string (testit.TestTester) ... MPI_COMM_WORLD size is 1
this, is, a, list, of, strings,
ok

----------------------------------------------------------------------
Ran 1 test in 0.147s

OK
```

## Installation

Do not pollute your python installation with this module.  Install it
in a personal directory. Here this directory is in the shell variable
`$prefix`. 

In order to complete the installation, `$prefix` must be in the
`PYTHONPATH` environment variable:
```
PYTHONPATH=$
export PYTHONPATH
```
or, if `PYTHONPATH` is not empty,

```
PYTHONPATH="${prefix}:$PYTHONPATH"
export PYTHONPATH
```
Install the module in `$prefix` with
```
python setup.py install --prefix="$prefix"
```
After installation, the `runme.py` script can be used
```
$prefix/bin/runme.py
```
which should produce output like this
```
MPI_COMM_WORLD size is 1
this, is, a, list, of, strings,

```

You should also be able to do something like

```
mpiexec -np 4 ~/Projects/GridPACK/gridpack-hadrec/bin/runme.py
```
and get something like

```
MPI_COMM_WORLD size is 4
this, is, a, list, of, strings,
MPI_COMM_WORLD size is 4
MPI_COMM_WORLD size is 4
this, is, a, list, of, strings, this, is, a, list, of, strings,

MPI_COMM_WORLD size is 4
this, is, a, list, of, strings,
```

