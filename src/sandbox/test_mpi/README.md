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
Or, 

```
python ./runme.py
```
This should produce output like

```
MPI_COMM_WORLD size is 1
```

## Installation

  
