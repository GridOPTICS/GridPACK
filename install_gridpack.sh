# This script installs GridPACK and python wrappter.GridPACK is built in src/build directory and installed in src/install directory.

# This script should be run from the top-level GridPACK directory.

# Flag for install GridPACK and GridPACK python wrapper
install_gridpack=true
install_gridpack_python=true

# Set your python executable here
python_exe=`which python`
if test -z ${python_exe}
then
    python_exe=`which python3`
fi

if test -z ${GRIDPACK_ROOT_DIR}
then
    export GRIDPACK_ROOT_DIR=${PWD}
    echo "GRIDPACK_ROOT_DIR = ${GRIDPACK_ROOT_DIR}"
fi

# Create directory for installing external packages
if test -z ${GP_EXT_DEPS}
then
     export GP_EXT_DEPS=${GRIDPACK_ROOT_DIR}/external-dependencies
fi

cd ${GRIDPACK_ROOT_DIR}

# Set environment variable GRIDPACK_BUILD_DIR and create a build directory
export GRIDPACK_BUILD_DIR=${GRIDPACK_ROOT_DIR}/src/build

GRIDPACK_INSTALL_DIR=${GRIDPACK_ROOT_DIR}/src/install
export GRIDPACK_DIR=${GRIDPACK_INSTALL_DIR}

if ${install_gridpack}
then

  rm -rf $GRIDPACK_BUILD_DIR
  mkdir $GRIDPACK_BUILD_DIR

  rm -rf ${GRIDPACK_INSTALL_DIR}

  cd ${GRIDPACK_BUILD_DIR}
    
  ## GridPACK installation
  echo "Building GridPACK"

#  git checkout develop

  rm -rf CMake*

  cmake_args="-D GA_DIR:STRING=${GP_EXT_DEPS}/ga-5.8/install_for_gridpack \
   -D BOOST_ROOT:STRING=${GP_EXT_DEPS}/boost_1_78_0/install_for_gridpack \    
   -D Boost_DIR:STRING=${GP_EXT_DEPS}/boost_1_78_0/install_for_gridpack/lib/cmake/Boost-1.78.0 \
   -D Boost_LIBRARIES:STRING=${GP_EXT_DEPS}/boost_1_78_0/install_for_gridpack/lib \   
   -D Boost_INCLUDE_DIRS:STRING=${GP_EXT_DEPS}/boost_1_78_0/install_for_gridpack/include \
   -D Boost_NO_BOOST_CMAKE:BOOL=TRUE \
   -D Boost_NO_SYSTEM_PATHS:BOOL=TRUE \
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

    rm -rf build
    
    ${python_exe} setup.py build

    rm -rf ${GRIDPACK_INSTALL_DIR}/lib/python
    mkdir ${GRIDPACK_INSTALL_DIR}/lib/python
    
    PYTHONPATH="${GRIDPACK_DIR}/lib/python:${PYTHONPATH}"
    export PYTHONPATH
    ${python_exe} setup.py install --home="$GRIDPACK_DIR"
    
fi

cd ${GRIDPACK_ROOT_DIR}

echo "Completed GridPACK installation"
