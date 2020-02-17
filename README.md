# GridPACK HADREC Application Python Wrapper

## Requirements

  * GridPACK >= 3.4
  * Python >= 2.7
  * pybind11 >= 2.4
  * Python `setuptools` package
  * Python `nose` package
  * CMake
  * C++ Compiler used to build GridPACK
  * MPI installation used to build GridPACK

### Pybind11

If pybind11 is not installed, or installed without CMake support, the
build will fail.  In that event, place the pybind11 code in the top
directory: 
```
git clone -b v2.4 https://github.com/pybind/pybind11.git
```
It is not needed after the module is built.
  
### GridPACK

GridPACK must be built and *installed* as *shared* libraries. This
requires that any GridPACK dependencies (e.g. PETSc, Global Arrays,
Boost) also be built as shared libraries. 

Set the `GRIDPACK_DIR` environment variable to indicate where GridPACK
was installed. 

## Build and Test

Build the GridPACK python wrapper modules with
```
GRIDPACK_DIR=/usr/local/gridpack
export GRIDPACK_DIR
python setup.py build
```
