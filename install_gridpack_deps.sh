
set -x

# This script installs all GridPACK dependencies.The dependencies are installed in external-dependencies directory.

# This script should be run from the top-level GridPACK directory.

# Flags for installing/not installing different dependency packages
install_boost=true
install_ga=true
install_petsc=true
install_shared=true

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

  boost_version="1.81.0"
  boost_us_version=`echo $boost_version | sed -e 's/\./_/g'`

  if "$install_shared"; then
      boostshared="shared"
  else
      boostshared="static"
  fi

  rm -rf boost*
    
  # Download and install Boost
  echo "Downloading Boost-$boost_version"

    # Download Boost
  wget https://boostorg.jfrog.io/artifactory/main/release/${boost_version}/source/boost_${boost_us_version}.tar.gz

  # Untar
  tar -xf boost_${boost_us_version}.tar.gz

  cd boost_${boost_us_version}

  # Install boost
  echo "Building Boost-${boost_version}"

  ./bootstrap.sh --prefix=install_for_gridpack --with-libraries=mpi,serialization,random,filesystem,system

  echo 'using mpi : mpicxx ; ' >> project-config.jam
  
  ./b2 -a -d+2 link="$boostshared" stage 

  echo "Installing Boost-$boost_version"
  ./b2 -a -d+2 link="$boostshared" install

  echo "Building and Installing Boost libraries complete"
fi

if ${install_ga}
then
  # Download, build, and install GA
  cd ${GP_EXT_DEPS}

  if "$install_shared"; then
      gashopts="--enable-shared=yes --enable-static=no"
  else
      gashopts="--enable-shared=no --enable-static=yes"
  fi

  ga_version="5.8.2"
  echo "Downloading GA-$ga_version"

  wget "https://github.com/GlobalArrays/ga/releases/download/v${ga_version}/ga-${ga_version}.tar.gz"

  tar -xf "ga-${ga_version}.tar.gz"

  cd ga-${ga_version}

  # Build GA
  echo "Building GA-${ga_version}"
  ./configure --with-mpi-ts --disable-f77 \
              --without-blas --without-lapack --without-scalapack \
              --enable-cxx --enable-i4 \
              --prefix=${PWD}/install_for_gridpack \
              CFLAGS=-Wno-implicit-function-declaration \
              $gashopts
  
  # Install GA
  echo "Installing GA-${ga_version}"
  make -j 10 install
  
  echo "GA-${ga_version} installation complete"
fi

if ${install_petsc}
then

  if "$install_shared"; then
      petscopts="--with-shared-libraries=1"
  else
      petscopts="--with-shared-libraries=0"
  fi

  # petsc_version="3.21.4"
  petsc_version="3.20.6"
  # petsc_version="3.19.4"
    
  # Install PETSc 
  cd ${GP_EXT_DEPS}

  rm -rf petsc*
  
  # Download
  echo "Downloading PETSc $petsc_version"

  git clone https://gitlab.com/petsc/petsc.git
    
  cd petsc

  # git checkout tags/v3.16.4 -b v3.16.4
  git checkout "tags/v$petsc_version" -b "v$petsc_version"

  export PETSC_DIR=${PWD}
  export PETSC_ARCH=build-dir
  rm -rf $PETSC_ARCH

  # Install PETSc
  echo "Installing PETSc $petsc_version"
    
  ./configure \
      --prefix=${PWD}/install_for_gridpack \
      --scalar-type=complex  \
      --with-fortran-bindings=0 \
      --download-superlu_dist \
      --download-metis \
      --download-parmetis \
      --download-suitesparse \
      --download-f2cblaslapack \
      --download-scalapack \
      --download-mumps \
      --download-cmake=0 \
      --with-sowing=0 \
      $petscopts

  # Build PETSc
  echo "Building PETSc $petsc_version"

  make 

  # Install PETSc
  echo "Installing PETSc $petsc_version"

  make install
  make check  
fi  

# update LD_LIBRARY_PATH so that boost,ga, and petsc are on it
export LD_LIBRARY_PATH=${GP_EXT_DEPS}/boost_${boost_us_version}/install_for_gridpack/lib:${GP_EXT_DEPS}/${ga_version}/install_for_gridpack/lib:${GP_EXT_DEPS}/petsc/install_for_gridpack/lib:${LD_LIBRARY_PATH}

# For Mac OS X
export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH

cd ${GRIDPACK_ROOT_DIR}

echo "Completed installing GridPACK dependencies in ${GP_EXT_DEPS}"
