# This script installs GridPACK and all its dependencies.The dependencies are installed in external-dependencies directory. GridPACK is built in src/build directory and installed in src/install directory.

# This script should be run from the top-level GridPACK directory.

# Set environment variable GRIDPACK_DIR
export GRIDPACK_DIR=${PWD}

# Create directory for installing external packages
rm -rf external-dependencies
mkdir external-dependencies

cd external-dependencies

export GP_EXT_DEPS=${PWD}

# Download and install Boost
echo "Downloading Boost-1.78.0"

# Download Boost
wget https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.gz

# Untar
tar -xvf boost_1_78_0.tar.gz

cd boost_1_78_0

# Install boost
echo "Building Boost-1.78.0"

./bootstrap.sh --prefix=install_for_gridpack --with-libraries=mpi,serialization,random,filesystem

cat project-config.jam >> `using mpi ;`

./b2 -a -d+2 link=shared stage

echo "Installing Boost-1.78.0"
./b2 -a -d+2 link=shared install

echo "Building and Installing Boost libraries complete"

# Download, build, and install GA
cd ..

echo "Downloading GA-5.8"

wget https://github.com/GlobalArrays/ga/releases/download/v5.8/ga-5.8.tar.gz

tar -xvf ga-5.8.tar.gz

cd ga-5.8

# Build GA
echo "Building GA-5.8"
./configure --with-mpi-ts --disable-f77 --without-blas --enable-cxx --enable-i4 --prefix=${PWD}/install_for_gridpack --enable-shared

# Install GA
echo "Installing GA-5.8"
make -j 10 install

echo "GA-5.8 installation complete"


# Install PETSc 3.16.4
cd ..

# Download
echo "Downloading PETSc 3.16.4"

git clone https://gitlab.com/petsc/petsc.githttps://gitlab.com/petsc/petsc.git

git checkout v3.16.4

cd petsc

export PETSC_DIR=${PWD}
export PETSC_ARCH=build-dir

# Install PETSc
echo "Installing PETSc 3.16.4"

./configure --download-superlu_dist --download-metis --download-parmetis --download-suitesparse --download-f2cblaslapack --download-cmake --prefix=${PWD}/install_for_gridpack --scalar-type=complex --with-shared-libraries=1 --download-f2cblaslapack=1 --download-cmake=1

# Build PETSc
echo "Building PETSc 3.16.4"

make

# Install PETSc
echo "Installing PETSc 3.16.4"

make install
make check

## GridPACK installation
echo "Building GridPACK develop branch"

git checkout develop

cd ${GRIDPACK_DIR}/src

rm -rf build
mkdir build

cd build

rm -rf CMake*                                                                               
cmake \                                                                                     
   -D GA_DIR:STRING=${GP_EXT_DEPS}/ga-5.8/install_for_gridpack \
   -D BOOST_ROOT:STRING=${GP_EXT_DEPS}/boost_1_78_0/install_for_gridpack \    
   -D Boost_DIR:STRING=${GP_EXT_DEPS}/boost_1_78_0/install_for_gridpack/lib/cmake/Boost-1.78.0 \   
   -D BOOST_LIBRARYDIR:STRING=${GP_EXT_DEPS}/boost_1_78_0/build/lib \   
   -D PETSC_DIR:PATH=${GP_EXT_DEPS}/petsc/install_for_gridpack \                  
   -D MPI_CXX_COMPILER:STRING='mpicxx' \     
   -D MPI_C_COMPILER:STRING='mpicc' \                                         
   -D MPIEXEC:STRING='mpiexec' \                                                            
   -D GRIDPACK_TEST_TIMEOUT:STRING=30 \                                                     
   -D CMAKE_INSTALL_PREFIX:PATH=${GRIDPACK_DIR}/src/install \             
   -D CMAKE_BUILD_TYPE:STRING=Debug \
   -D BUILD_SHARED_LIBS=YES \
    ..   

echo "Installing GridPACK develop branch"
make -j 10 install
