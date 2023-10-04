#! /bin/bash

# installs GridPACK and its python wrapper
# gridPACK is built in src/build and installed to src/install
# run this from the top level GridPACK directory

set -xeuo pipefail

function install_gridpack {
  echo "--- GridPACK ---"

  # args
  local gridpack_deps_dir=$1
  local gridpack_build_dir=$2
  local gridpack_install_dir=$3
  : "${gridpack_deps_dir:?}" "${gridpack_build_dir:?}" "${gridpack_install_dir:?}"

  # remove existing build and install dir
  rm -rf "$gridpack_build_dir"
  rm -rf "$gridpack_install_dir"

  # create the build dir
  mkdir -p "$gridpack_build_dir"

  pushd "$gridpack_build_dir" || exit

  # todo? git checkout develop

  # remove existing cmake output
  rm -rf CMake*

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
    -D GRIDPACK_TEST_TIMEOUT:STRING=30 \
    -D CMAKE_INSTALL_PREFIX:PATH="${gridpack_install_dir}" \
    -D CMAKE_BUILD_TYPE:STRING=Debug \
    -D BUILD_SHARED_LIBS=YES \
    -D Boost_NO_SYSTEM_PATHS:BOOL=TRUE \
    -D Boost_NO_BOOST_CMAKE:BOOL=TRUE \
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
  local gridpack_build_dir=$1
  local gridpack_install_dir=$2
  : "${gridpack_build_dir:?}" "${gridpack_install_dir:?}"

  # update submodules
  echo "Updating GridPACK submodules"
  git submodule update --init

  pushd python || exit

  os_id=$(
    source /etc/os-release
    echo "$ID"
  )
  if [[ $os_id == "rhel" ]] || [[ $os_id == "centos" ]]; then
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

echo "Installing GridPACK"
date

build_dir=${PWD}/src/build
install_dir=${PWD}/src/install

install_gridpack "${GP_EXT_DEPS:?}" "$build_dir" "$install_dir"
install_gridpack_python "$build_dir" "$install_dir"

echo "Completed GridPACK installation"
