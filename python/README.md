# GridPACK Python Wrapper

This provides Python wrappers for the
[GridPACK](https://github.com/GridOPTICS/GridPACK/) library: power flow,
dynamic simulation, state estimation, HADREC, EMT and contingency analysis,
plus a `gridpack` command-line interface for all six.

## Quick start (Docker)

`gridpack` is a compiled extension linked against GridPACK's shared
libraries, so it needs GridPACK, PETSc, Global Arrays, Boost and MPI already
built.  The published image has all of that, so this adds only the Python
layer and takes minutes instead of hours:

```
docker build -f ../Dockerfile.pygridpack -t pygridpack .
docker run --rm -it -v "$PWD:/app/workspace" pygridpack
```

Inside the container, no environment variables are needed.  As a CLI:

```
gridpack powerflow input_14.xml
mpiexec -np 4 gridpack powerflow input_14.xml
```

Or as a library:

```python
from gridpack import Session, PowerFlow

with Session() as s:                       # one Session per process
    result = PowerFlow(s, "input_14.xml").solve()
    print(result.converged)
```

Run the test suite against the container's own input data:

```
docker run --rm pygridpack pytest /app/python/tests -q
```

To build from source instead, read on.

## Requirements

  * GridPACK >= 3.4, built and *installed* as shared libraries
  * Python >= 3.9, < 3.14
  * pybind11 >= 2.4 (included as a submodule, see below)
  * `pip` and `setuptools` >= 61
  * `mpi4py`, built against the same MPI as GridPACK
  * The CMake, C++ compiler and MPI installation used to build GridPACK
  * `pytest`, to run the test suite

### Pybind11

The source for `pybind11` is needed to build the wrapper, so it is included
as a submodule.  Make sure it is present before building:
```
git submodule update --init
```

### GridPACK

GridPACK must be built and *installed* as *shared* libraries. This
requires that any GridPACK dependencies (e.g. PETSc, Global Arrays,
Boost) also be built as shared libraries.

## Build and install from source

### 1. Point at the GridPACK installation

`GRIDPACK_DIR` is how the build finds GridPACK:
```
export GRIDPACK_DIR=/path/to/gridpack/install
```

### 2. Install into a virtual environment

A virtual environment is the recommended target: the module lands somewhere
Python already looks, so no `PYTHONPATH` is needed.
```
python -m venv ~/.venvs/pygridpack
source ~/.venvs/pygridpack/bin/activate
pip install mpi4py
pip install --no-build-isolation --no-deps -e python/
```
`--no-build-isolation` matters.  `mpi4py` is a build requirement, and an
isolated build would fetch and compile its own copy against whichever MPI it
happens to find rather than the one GridPACK was built with.  Drop `-e` for a
non-editable install, or install `python/[test]` instead to get the test
dependencies too.

### 3. Check

```
python -c 'import gridpack; print(gridpack.__all__)'
gridpack --version
```

### Installing alongside GridPACK instead

To put the module in the same prefix as the rest of GridPACK:
```
pip install --no-build-isolation --no-deps --upgrade --prefix=$GRIDPACK_DIR python/
```
`PYTHONPATH` then has to include the site-packages directory under that
prefix, which depends on the Python version:
```
export PYTHONPATH=$GRIDPACK_DIR/lib/python3.12/site-packages
```

## Test

The `test` extra pulls in pytest, along with the `pandas` and `numpy` the
result helpers use:
```
pip install --no-build-isolation -e "python/[test]"
python -m pytest python/tests -q
```
Tests needing a working GridPACK build are marked `integration`, and those
needing `mpiexec` are marked `mpi`, so the pure-Python surface can be checked
on its own:
```
python -m pytest python/tests -q -m "not integration"
```

## Run examples

`python/src/example` holds scripts that drive the interface directly:
```
cd python/src/example
python 39bus_test_example.py
python 39bus_test_example_dsf.py
python 39bus_scatterload_steptest_new_itr.py
mpiexec -np 2 python 39bus_scatterload_steptest_new_itr.py
python 39bus_scatterload_steptest_new_itr_dsf.py
mpiexec -np 2 python 39bus_scatterload_steptest_new_itr_dsf.py
python 39bus_scatterload_steptest_new_itr_compensateY.py
python 39bus_test_pfdata.py
```

Or use the `gridpack` CLI, which covers all six applications:
```
gridpack --help
gridpack powerflow case.xml
mpiexec -np 4 gridpack powerflow case.xml
gridpack ca case.xml --vlimits 0.9 1.1     # needs <contingencyList> in the XML
```
GridPACK resolves `<networkConfiguration>` relative to the working directory,
so the network file has to sit beside the XML.  A missing one aborts with
`wrong dimension specified` rather than a file-not-found error, so check that
first.  The stock inputs keep the XML and the `.raw` in separate directories
and need staging into one:
```
export GP=/path/to/GridPACK          # the source tree
mkdir /tmp/pf14 && cd /tmp/pf14
cp $GP/src/applications/data_sets/input/powerflow/input_14.xml \
   $GP/src/applications/data_sets/raw/IEEE14.raw .
gridpack powerflow input_14.xml
mpiexec -np 4 gridpack powerflow input_14.xml
```
The CLI exits 0 on success, 2 on bad usage or an unreadable config, 3 when the
analysis did not converge, and 1 on anything unexpected, so it can be driven
from a script.
