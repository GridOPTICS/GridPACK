#! /bin/bash

# installs GridPACK and its python wrapper
# gridPACK is built in src/build and installed to src/install
# run this from the top level GridPACK directory

# bash options:
# - xtrace: print each command before executing it
# - errexit: exit on error
# - nounset: treat unset variables as errors
# - pipefail: treat whole pipeline as errored if any commands within error
# https://www.gnu.org/software/bash/manual/html_node/The-Set-Builtin.html
set -o xtrace -o errexit -o nounset -o pipefail

function install_gridpack {
  echo "--- GridPACK ---"

  # args
  local gridpack_deps_dir=${1:?}
  local gridpack_build_dir=${2:?}
  local gridpack_install_dir=${3:?}

  # remove existing build and install dir
  rm -rf "$gridpack_build_dir"
  rm -rf "$gridpack_install_dir"

  # create the build dir
  mkdir -p "$gridpack_build_dir"

  pushd "$gridpack_build_dir" || exit

  # todo? git checkout develop

  # remove existing cmake output
  rm -rf CMake*

  # load mpi module for RHEL
  load_mpi_module

  # generate make files
  echo "Generating GridPACK make files"
  cmake \
    -D GA_DIR:STRING="${gridpack_deps_dir}/ga/install_for_gridpack" \
    -D BOOST_ROOT:STRING="${gridpack_deps_dir}/boost/install_for_gridpack" \
    -D Boost_DIR:STRING="${gridpack_deps_dir}/boost/install_for_gridpack/lib/cmake/Boost" \
    -D Boost_LIBRARIES:STRING="${gridpack_deps_dir}/boost/install_for_gridpack/lib" \
    -D Boost_INCLUDE_DIRS:STRING="${gridpack_deps_dir}/boost/install_for_gridpack/include" \
    -D PETSC_DIR:PATH="${gridpack_deps_dir}/petsc/install_for_gridpack" \
    -D MPI_CXX_COMPILER:STRING='mpicxx' \
    -D MPI_C_COMPILER:STRING='mpicc' \
    -D MPIEXEC:STRING='mpiexec' \
    -D GRIDPACK_TEST_TIMEOUT:STRING=120 \
    -D CMAKE_INSTALL_PREFIX:PATH="${gridpack_install_dir}" \
    -D CMAKE_BUILD_TYPE:STRING=Debug \
    -D BUILD_SHARED_LIBS=YES \
    -D Boost_NO_SYSTEM_PATHS:BOOL=TRUE \
    -D Boost_NO_BOOST_CMAKE:BOOL=TRUE \
    -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
    ..

  # install
  echo "Installing GridPACK"
  make -j "${MAKE_JOBS:-$(nproc)}" install

  popd || exit

  echo "GridPACK installation complete"
}

function install_gridpack_python {
  echo "--- GridPACK python wrapper ---"

  # args
  local gridpack_build_dir=${1:?}
  local gridpack_install_dir=${2:?}

  # update submodules
  echo "Updating GridPACK submodules"
  git submodule update --init

  pushd python || exit

  # export RHEL_OPENMPI_HACK for RHEL
  distribution=$(get_base_distro)
  if test "$distribution" = "rhel"; then
    export RHEL_OPENMPI_HACK=yes
  fi

  # set python executable path
  python_exe=$(which python || which python3)

  # remove existing build dir
  rm -rf build

  # export GRIDPACK_DIR
  export GRIDPACK_DIR="${gridpack_install_dir}"

  # build
  echo "Building GridPACK python wrapper"
  ${python_exe} setup.py build

  # remove existing install dir
  local py_lib="${gridpack_install_dir}/lib/python"
  rm -rf "${py_lib}"
  mkdir -p "${py_lib}"

  # add lib to python path
  export PYTHONPATH="${py_lib}${PYTHONPATH:+:$PYTHONPATH}"

  # install
  echo "Installing GridPACK python wrapper"
  ${python_exe} setup.py install --home="$gridpack_install_dir"

  popd || exit

  echo "GridPACK python wrapper installation complete"
}

# get the parent directory of this script
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

# load the `load_mpi_module`, `get_base_distro` functions
source "$script_dir/install_package_deps_lib.sh"

echo "Installing GridPACK"
date

build_dir=${PWD}/src/build
install_dir=${PWD}/src/install

gp_ext_deps=${GP_EXT_DEPS:-/gridpack-dependencies}

install_gridpack "$gp_ext_deps" "$build_dir" "$install_dir"
install_gridpack_python "$build_dir" "$install_dir"

echo "Completed GridPACK installation"
