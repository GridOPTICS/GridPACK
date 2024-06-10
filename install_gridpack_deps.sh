# This script installs all GridPACK dependencies.The dependencies are installed in external-dependencies directory.

# This script should be run from the top-level GridPACK directory.

# Flags for installing/not installing different dependency packages
install_boost=true
install_ga=true
install_petsc=true

echo $(date)

if test -z ${GRIDPACK_ROOT_DIR}
then
    export GRIDPACK_ROOT_DIR=${PWD}
    echo "GRIDPACK_ROOT_DIR = ${GRIDPACK_ROOT_DIR}"
fi

cd ${GRIDPACK_ROOT_DIR}

# Create directory for installing external packages
if test -z ${GP_EXT_DEPS}
then
     export GP_EXT_DEPS=${GRIDPACK_ROOT_DIR}/external-dependencies
     rm -rf ${GP_EXT_DEPS}
     mkdir ${GP_EXT_DEPS}
     echo "GRIDPACK_EXT_DEPENDENCIES_DIR=${GP_EXT_DEPS}"
else
    if test -d ${GP_EXT_DEPS}
    then
	echo "GRIDPACK_EXT_DEPENDENCIES_DIR=${GP_EXT_DEPS}"
    else
	mkdir ${GP_EXT_DEPS}
    fi		
fi

cd ${GP_EXT_DEPS}

if ${install_boost}
then

  rm -rf boost*
    
  # Download and install Boost
  echo "Downloading Boost-1.81.0"

  # Download Boost
  wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz

  # Untar
  tar -xf boost_1_81_0.tar.gz

  cd boost_1_81_0

  # Install boost
  echo "Building Boost-1.81.0"

  ./bootstrap.sh --prefix=install_for_gridpack --with-libraries=mpi,serialization,random,filesystem,system

  echo 'using mpi ;' >> project-config.jam

  ./b2 -a -d+2 link=static stage

  echo "Installing Boost-1.78.0"
  ./b2 -a -d+2 link=static install

  echo "Building and Installing Boost libraries complete"
fi

if ${install_ga}
then
  # Download, build, and install GA
  cd ${GP_EXT_DEPS}

  rm -rf ga*

  echo "Downloading GA-5.8"

  wget https://github.com/GlobalArrays/ga/releases/download/v5.8/ga-5.8.tar.gz

  tar -xf ga-5.8.tar.gz

  cd ga-5.8

  # Build GA
  echo "Building GA-5.8"
  ./configure --with-mpi-ts --disable-f77 --without-blas --enable-cxx --enable-i4 --prefix=${PWD}/install_for_gridpack
  
  # Install GA
  echo "Installing GA-5.8"
  make -j 10 install
  
  echo "GA-5.8 installation complete"
fi

if ${install_petsc}
then
    
  # Install PETSc 3.20.1
  cd ${GP_EXT_DEPS}

  rm -rf petsc*
  
  # Download
  echo "Downloading PETSc"

  git clone https://gitlab.com/petsc/petsc.git
    
  cd petsc

  git checkout tags/v3.20.1 -b v3.20.1

  export PETSC_DIR=${PWD}
  export PETSC_ARCH=build-dir

  # Install PETSc
  echo "Installing PETSc 3.20.1"
    
  ./configure --download-superlu_dist --download-metis --download-parmetis --with-debugging=0 --prefix=${PWD}/install_for_gridpack

  # Build PETSc
  echo "Building PETSc 3.20.1"

  make 

  # Install PETSc
  echo "Installing PETSc 3.20.1"

  make install
  make check  
fi  

# update LD_LIBRARY_PATH so that boost,ga, and petsc are on it
export LD_LIBRARY_PATH=${GP_EXT_DEPS}/boost_1_81_0/install_for_gridpack/lib:${GP_EXT_DEPS}/ga-5.8/install_for_gridpack/lib:${GP_EXT_DEPS}/petsc/install_for_gridpack/lib:${LD_LIBRARY_PATH}

cd ${GRIDPACK_ROOT_DIR}

echo "Completed installing GridPACK dependencies in ${GP_EXT_DEPS}"
