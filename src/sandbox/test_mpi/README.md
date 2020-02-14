# Simple PyBind11 Wrapper

## Requirements

  * Python >= 2.7
  * pybind11 >= 2.4
  * Python `setuptools` package
  * Python `nose` package
  * CMake
  * C++ Compiler
  * MPI 
  
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
