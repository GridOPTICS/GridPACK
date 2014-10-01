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
    CFLAGS="-pthread"
    export CFLAGS
    CXXFLAGS="-pthread"
    export CXXFLAGS
    FC="$prefix/bin/gfortran"
    export FC
    FCFLAGS="-pthread"
    export FCFLAGS

    cmake -Wdev --debug-trycompile \
        -D CMAKE_Fortran_FLAGS:STRING="-pthread" \
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

    CC=/share/apps/gcc/4.6.2/bin/gcc
    export CC
    CXX=/share/apps/gcc/4.6.2/bin/g++
    export CXX
    FC=/share/apps/gcc/4.6.2/bin/gfortran
    export FC

    prefix="/pic/projects/gridpack/olympus-gnu-4.6"
    cmake -Wdev --debug-trycompile \
        -D GRIDPACK_DIR:PATH="$prefix/gridpack" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack" \
        $common_flags ..

else

    echo "Unknown host: $host"
    exit 2
    
fi

