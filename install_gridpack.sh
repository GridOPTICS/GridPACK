# This script installs GridPACK and all its dependencies.The dependencies are installed in src/build/external-dependencies directory. GridPACK is built in src/build directory and installed in src/install directory.

# This script should be run from the top-level GridPACK directory.

# Flags for installing/not installing different dependency packages
install_boost=true
install_ga=true
install_petsc=true

# Flag for install GridPACK and GridPACK python wrapper
install_gridpack=true
install_gridpack_python=true

# Set your python executable here
python_exe=`which python`
if[ -z ${python_exe}]
then
    python_exe=`which python3`
fi

# Install log file
install_log_file=${PWD}/install.log

echo 'GRIDPACK installation log file' > ${install_log_file}
echo $(date) >> install_log_file

export GRIDPACK_ROOT_DIR=${PWD}
# Set environment variable GRIDPACK_BUILD_DIR and create a build directory
export GRIDPACK_BUILD_DIR=${GRIDPACK_ROOT_DIR}/src/build

rm -rf $GRIDPACK_BUILD_DIR
mkdir $GRIDPACK_BUILD_DIR
cd $GRIDPACK_BUILD_DIR

# Create directory for installing external packages
mkdir external-dependencies

export GP_EXT_DEPS=${GRIDPACK_BUILD_DIR}/external-dependencies

cd external-dependencies

if ${install_boost}
then

  rm -rf boost*
    
  # Download and install Boost
  echo "Downloading Boost-1.78.0"

  # Download Boost
  wget https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.gz &>> ${install_log_file}

  # Untar
  tar -xf boost_1_78_0.tar.gz

  cd boost_1_78_0

  # Install boost
  echo "Building Boost-1.78.0"

  ./bootstrap.sh --prefix=install_for_gridpack --with-libraries=mpi,serialization,random,filesystem,system &>> ${install_log_file}

  echo 'using mpi ;' >> project-config.jam

  ./b2 -a -d+2 link=shared stage &>> ${install_log_file}

  echo "Installing Boost-1.78.0"
  ./b2 -a -d+2 link=shared install &>> ${install_log_file}

  echo "Building and Installing Boost libraries complete"
fi

if ${install_ga}
then
  # Download, build, and install GA
  cd ${GP_EXT_DEPS}

  echo "Downloading GA-5.8"

  wget https://github.com/GlobalArrays/ga/releases/download/v5.8/ga-5.8.tar.gz  &>> ${install_log_file}

  tar -xf ga-5.8.tar.gz  &>> ${install_log_file}

  cd ga-5.8

  # Build GA
  echo "Building GA-5.8"
  ./configure --with-mpi-ts --disable-f77 --without-blas --enable-cxx --enable-i4 --prefix=${PWD}/install_for_gridpack --enable-shared &>> ${install_log_file}
  
  # Install GA
  echo "Installing GA-5.8"
  make -j 10 install &>> ${install_log_file}
  
  echo "GA-5.8 installation complete"
fi

if ${install_petsc}
then
    
  # Install PETSc 3.16.4
  cd ${GP_EXT_DEPS}

  # Download
  echo "Downloading PETSc 3.16.4"

  git clone https://gitlab.com/petsc/petsc.git &>> ${install_log_file}
    
  cd petsc

  git checkout tags/v3.16.4 -b v3.16.4  &>> ${install_log_file}

  export PETSC_DIR=${PWD}
  export PETSC_ARCH=build-dir

  # Install PETSc
  echo "Installing PETSc 3.16.4"
    
  ./configure --download-superlu_dist --download-metis --download-parmetis --download-suitesparse --download-f2cblaslapack --download-cmake --prefix=${PWD}/install_for_gridpack --scalar-type=complex --with-shared-libraries=1 --download-f2cblaslapack=1  &>> ${install_log_file}

  # Build PETSc
  echo "Building PETSc 3.16.4"

  make  &>> ${install_log_file}

  # Install PETSc
  echo "Installing PETSc 3.16.4"

  make install  &>> ${install_log_file}
  make check  &>> ${install_log_file}
fi  

GRIDPACK_INSTALL_DIR=${GRIDPACK_ROOT_DIR}/src/install
export GRIDPACK_DIR=${GRIDPACK_INSTALL_DIR}

if ${install_gridpack}
then

  cd ${GRIDPACK_BUILD_DIR}
    
  ## GridPACK installation
  echo "Building GridPACK develop branch"

  git checkout develop

  rm -rf CMake*

  cmake_args="-D GA_DIR:STRING=${GP_EXT_DEPS}/ga-5.8/install_for_gridpack \
   -D BOOST_ROOT:STRING=${GP_EXT_DEPS}/boost_1_78_0/install_for_gridpack \    
   -D Boost_DIR:STRING=${GP_EXT_DEPS}/boost_1_78_0/install_for_gridpack/lib/cmake/Boost-1.78.0 \   
   -D BOOST_LIBRARIES:STRING=${GP_EXT_DEPS}/boost_1_78_0/build/lib \   
   -D PETSC_DIR:PATH=${GP_EXT_DEPS}/petsc/install_for_gridpack \                  
   -D MPI_CXX_COMPILER:STRING='mpicxx' \     
   -D MPI_C_COMPILER:STRING='mpicc' \                                         
   -D MPIEXEC:STRING='mpiexec' \                                                            
   -D GRIDPACK_TEST_TIMEOUT:STRING=30 \                                                     
   -D CMAKE_INSTALL_PREFIX:PATH=${GRIDPACK_INSTALL_DIR} \             
   -D CMAKE_BUILD_TYPE:STRING=Debug \
   -D BUILD_SHARED_LIBS=YES \
   -D Boost_NO_SYSTEM_PATHS:BOOL=TRUE \
    ..   "
 
  cmake ${cmake_args}                                                                                     

   echo "Installing GridPACK develop branch"
   make -j 10 install
fi

if ${install_gridpack_python}
then
    echo "Installing GridPACK python wrapper"
    
    cd ${GRIDPACK_ROOT_DIR}
    
    git submodule update --init

    echo ${GRIDPACK_DIR}
    cd python

    export RHEL_OPENMPI_HACK=yes

    ${python_exe} setup.py build

    rm -rf ${GRIDPACK_INSTALL_DIR}/lib/python
    mkdir ${GRIDPACK_INSTALL_DIR}/lib/python
    
    PYTHONPATH="${GRIDPACK_DIR}/lib/python:${PYTHONPATH}"
    export PYTHONPATH
    ${python_exe} setup.py install --home="$GRIDPACK_DIR"
    
fi

# update LD_LIBRARY_PATH so that boost,ga, and petsc are on it
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GP_EXT_DEPS}/boost_1_78_0/install_for_gridpack/lib:${GP_EXT_DEPS}/ga-5.8/install_for_gridpack/lib:${GP_EXT_DEPS}/petsc/install_for_gridpack/lib
cd ${GRIDPACK_ROOT_DIR}

echo "Completed GridPACK installation"
