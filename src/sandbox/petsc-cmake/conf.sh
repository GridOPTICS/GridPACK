#! /bin/sh

host=`uname -n`

rm -f CMakeCache.txt 

if [ $host == "flophouse" ]; then

    prefix="/net/flophouse/files0/perksoft/linux64"
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

    cmake -Wno-dev \
        -D Boost_DIR:STRING="$prefix" \
        -D PETSC_DIR:STRING="/net/flophouse/files0/perksoft/petsc-3.4.3" \
        -D PETSC_ARCH:STRING='arch-linux2-complex-opt' \
        -D MPI_CXX_COMPILER:STRING="$prefix/bin/mpicxx" \
        -D MPI_C_COMPILER:STRING="$prefix/bin/mpicc" \
        -D MPIEXEC:STRING="$prefix/bin/mpiexec" \
        -D CMAKE_BUILD_TYPE:STRING="Debug" \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        ..
    
elif [ $host == "pe10900" ]; then

    cmake -Wno-dev \
        -D Boost_DIR:STRING='/opt/local' \
        -D PETSC_DIR:STRING='/net/flophouse/files0/perksoft/petsc-3.3-p3' \
        -D PETSC_ARCH:STRING='arch-darwin-cxx-opt' \
        -D MPI_CXX_COMPILER:STRING='openmpicxx' \
        -D MPI_C_COMPILER:STRING='openmpicc' \
        -D MPIEXEC:STRING='openmpiexec' \
        -D CMAKE_BUILD_TYPE:STRING="Debug" \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        ..

else

    echo "Unknown host: $host"
    exit 2
    
fi

