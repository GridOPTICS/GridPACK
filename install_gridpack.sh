# This script installs GridPACK and python wrappter.GridPACK is built in src/build directory and installed in src/install directory.

set -x

# This script should be run from the top-level GridPACK directory.

# Flag for install GridPACK and GridPACK python wrapper
install_gridpack=true
install_gridpack_shared=false
install_gridpack_python=true

# These must match install_gridpack_deps.sh
boost_version="1.81.0"
ga_version="5.8.2"

if "$install_gridpack_python"; then
    install_gridpack_shared=true
    # Set your python executable here
    python_exe=`which python`
    if test -z ${python_exe}
    then
        python_exe=`which python3`
    fi
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

boost_us_version=`echo $boost_version | sed -e 's/\./_/g'`
boost_dir="${GP_EXT_DEPS}/boost_${boost_us_version}/install_for_gridpack"

petsc_dir="${GP_EXT_DEPS}/petsc/install_for_gridpack"

ga_dir="${GP_EXT_DEPS}/ga-${ga_version}/install_for_gridpack"

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

  cmake_args="-D GA_DIR:STRING=${ga_dir}
   -D BOOST_ROOT:STRING=$boost_dir 
   -D Boost_DIR:STRING=$boost_dir 
   -D Boost_LIBRARIES:STRING=${boost_dir}/lib 
   -D Boost_INCLUDE_DIRS:STRING=${boost_dir}/include 
   -D Boost_NO_BOOST_CMAKE:BOOL=TRUE 
   -D Boost_NO_SYSTEM_PATHS:BOOL=TRUE 
   -D PETSC_DIR:PATH=${petsc_dir} 
   -D MPI_CXX_COMPILER:STRING='mpicxx' 
   -D MPI_C_COMPILER:STRING='mpicc' 
   -D MPIEXEC:STRING='mpiexec' 
   -D MPIEXEC_MAX_NUMPROCS:STRING=2
   -D GRIDPACK_TEST_TIMEOUT:STRING=120 
   -D CMAKE_INSTALL_PREFIX:PATH=${GRIDPACK_INSTALL_DIR} 
   -D CMAKE_BUILD_TYPE:STRING=Debug 
   -D BUILD_SHARED_LIBS=$install_gridpack_shared 
    ..   "
 
  cmake -Wdev ${cmake_args}                                                                                     

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

    # This is necessary on RHEL/CentOS systems
    # export RHEL_OPENMPI_HACK=yes

    rm -rf build dist gridpack_hadrec.egg-info/
    
    ${python_exe} -m pip install --no-deps --upgrade --prefix=$GRIDPACK_INSTALL_DIR .

    # Construct the installed package path and check it. The actual
    # path with in the prefix seems to be kind of random, so
    # add a couple of options to PYTHONPATH

    pyvnum=`${python_exe} -V | sed -e 's/Python  *\([23]\)\.\([0-9][0-9]*\)\..*/\1.\2/' `
    pydir="${GRIDPACK_DIR}/lib/python${pyvnum}/site-packages:${GRIDPACK_DIR}/local/lib/python${pyvnum}/dist-packages"
    
    if [ -z "$PYTHONPATH" ]; then
        PYTHONPATH="${pydir}:${PYTHONPATH}"
    else
        PYTHONPATH="${pydir}"
    fi
    export PYTHONPATH

    # A quick check to see if the Python module is available
    ${python_exe} -c 'import gridpack'
    
fi

cd ${GRIDPACK_ROOT_DIR}

echo "Completed GridPACK installation"

set +x
