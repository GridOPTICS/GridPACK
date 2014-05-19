#! /bin/sh

host=`uname -n`

rm -f CMakeCache.txt 

options="-Wdev --debug-trycompile"

# useful build types: Debug, Release, RelWithDebInfo
common_flags="\
        -D CMAKE_BUILD_TYPE:STRING=Debug \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
"

if [ $host == "flophouse" ]; then

    CC=/usr/bin/gcc44
    export CC
    CXX=/usr/bin/g++44
    export CXX
    CFLAGS="-pthread"
    export CFLAGS
    CXXFLAGS="-pthread"
    export CXXFLAGS
    FC=/usr/bin/gfortran44
    export FC

    prefix="/net/flophouse/files0/perksoft/linux64/openmpi44"
    cmake -Wdev --debug-trycompile \
        -D CLASS_WORD:STRING="TYPE" \
        -D GRIDPACK_DIR:PATH="$prefix/gridpack" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack" \
        $common_flags ..
    
elif [ $host == "pe10900" ]; then
    prefix="/net/flophouse/files0/perksoft/macosx"

    # avoid using the system compilers and MPI wrappers -- use MacPorts

    CC=/opt/local/bin/gcc
    export CC
    CXX=/opt/local/bin/g++
    export CXX
    FC=/opt/local/bin/gfortran
    export FC

    cmake $options \
        -D GRIDPACK_DIR:PATH="$prefix/gridpack" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack" \
        $common_flags ..

elif [ $host == "olympus.local" ]; then

    prefix="/pic/projects/gridpack/software"
    cmake $options \
        -D GA_DIR:STRING="/pic/projects/gridpack/ga-5-2" \
        -D GA_EXTRA_LIBS:STRING="-libverbs" \
	-D BOOST_ROOT:STRING="$prefix" \
	-D PETSC_DIR:STRING="$prefix/petsc-3.4.0" \
	-D PETSC_ARCH:STRING='olympus-openmpi-gnu-cxx-complex-opt' \
	-D MPI_CXX_COMPILER:STRING='mpicxx' \
	-D MPI_C_COMPILER:STRING='mpicc' \
	-D MPIEXEC:STRING='mpiexec' \
	$common_flags ..

else

    echo "Unknown host: $host"
    exit 2
    
fi

