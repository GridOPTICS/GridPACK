#! /bin/bash

# get the parent directory of this script
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

petsc_version="${PETSC_VERSION:-3.16.4}"

echo "--- Installing PETSc ${petsc_version} ---"

# clone
echo "Cloning PETSc repository"
git clone https://gitlab.com/petsc/petsc.git

pushd petsc || exit

git checkout "tags/v${petsc_version}" -b "v${petsc_version}"

export PETSC_DIR=${PWD}
export PETSC_ARCH=build-dir

# load mpi module for RHEL
source "$script_dir/install_package_deps_lib.sh"
load_mpi_module

# install
echo "Configuring PETSc"
./configure \
  --download-mumps \
  --download-scalapack \
  --download-metis \
  --download-parmetis \
  --download-suitesparse \
  --download-f2cblaslapack \
  --prefix="${PWD}"/install_for_gridpack \
  --scalar-type=complex \
  --with-shared-libraries=1

# build
echo "Building PETSc"
make

# install
echo "Installing PETSc"
make install

# check
echo "Checking PETSc"
make check

popd || exit

echo "PETSc installation complete"
