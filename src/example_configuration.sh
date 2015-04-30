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


    prefix="/net/flophouse/files0/perksoft/linux64/openmpi48"
    PATH="${prefix}/bin:${PATH}"
    export PATH

    CC="$prefix/bin/gcc"
    export CC
    CXX="$prefix/bin/g++"
    export CXX
    CFLAGS="-pthread -Wall"
    export CFLAGS
    CXXFLAGS="-pthread -Wall"
    export CXXFLAGS

    cmake -Wdev --debug-trycompile \
        -D GA_DIR:STRING="$prefix" \
        -D BOOST_ROOT:STRING="$prefix" \
        -D PETSC_DIR:STRING="/net/flophouse/files0/perksoft/petsc-3.5.2" \
        -D PETSC_ARCH:STRING='linux-gnu48-complex-opt' \
        -D MPI_CXX_COMPILER:STRING="$prefix/bin/mpicxx" \
        -D MPI_C_COMPILER:STRING="$prefix/bin/mpicc" \
        -D MPIEXEC:STRING="$prefix/bin/mpiexec" \
        -D MPIEXEC_MAX_NUMPROCS:STRING="4" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack" \
        $common_flags ..
    
elif [ $host == "pe10900" ]; then
    prefix="/net/flophouse/files0/perksoft/macosx"

    # avoid using the system compilers and MPI wrappers -- use MacPorts

    CC=/opt/local/bin/gcc
    CFLAGS="-Wreturn-type"
    export CC CFLAGS
    CXX=/opt/local/bin/g++
    CXXFLAGS="-Wreturn-type"
    export CXX CXXFLAGS

    cmake $options \
        -D GA_DIR:STRING="$prefix" \
        -D GA_EXTRA_LIBS:STRING="-lblas" \
        -D BOOST_ROOT:STRING='/opt/local' \
        -D PETSC_DIR:STRING="$prefix/../petsc-3.4.3" \
        -D PETSC_ARCH:STRING='arch-macosx-complex-opt' \
        -D MPI_CXX_COMPILER:STRING='/opt/local/bin/mpicxx' \
        -D MPI_C_COMPILER:STRING='/opt/local/bin/mpicc' \
        -D MPIEXEC:STRING='/opt/local/bin/mpiexec' \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING="10" \
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

elif [ $host == "gridpack1" ]; then

    CC=/usr/bin/gcc
    export CC
    CXX=/usr/bin/g++
    export CXX
    CFLAGS="-pthread"
    export CFLAGS
    CXXFLAGS="-pthread"
    export CXXFLAGS

    prefix="/net/flophouse/files0/perksoft/linux64/openmpi44"
    cmake -Wdev --debug-trycompile \
        -D GA_DIR:STRING="/usr/local/ga-5-3" \
        -D BOOST_ROOT:STRING="/home/d3m998/boost_1_53_0" \
        -D PETSC_DIR:STRING="/home/d3m998/Software/petsc-3.4.5" \
        -D PETSC_ARCH:STRING='arch-linux2-complex-opt' \
        -D MPI_CXX_COMPILER:STRING="/usr/local/openmpi/bin/mpicxx" \
        -D MPI_C_COMPILER:STRING="/usr/local/openmpi/bin/mpicc" \
        -D MPIEXEC:STRING="/usr/local/openmpi/bin/mpiexec" \
        -D CMAKE_INSTALL_PREFIX:PATH="/home/d3m998/gridpack-Github/install" \
        $common_flags ..

else

    echo "Unknown host: $host"
    exit 2
    
fi

