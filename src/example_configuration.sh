#! /bin/sh


# -------------------------------------------------------------
# handle command line options
# -------------------------------------------------------------
usage="$0 [-d|-r] [name]"

set -- `getopt d $*`
if [ $? != 0 ]; then
    echo $usage >&2
    exit 2
fi

build="RelWithDebInfo"
for o in $*; do
    case $o in
        -d)
            build="Debug"
            shift
            ;;
        -r)
            build="Release"
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "$0: error: $o: unknown option" >&2
            echo $usage >&2
            exit 2
    esac
done

if [ $# -gt 0 ]; then
    host="$1"
else
    host=`uname -n`
fi

rm -rf CMakeCache.txt CMakeFiles

options="-Wdev --debug-trycompile"

# useful build types: Debug, Release, RelWithDebInfo
common_flags="\
        -D CMAKE_BUILD_TYPE:STRING=$build \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
"

if [ $host == "flophouse" ]; then

    # RHEL 5 with GNU 4.8 compilers built from scratch

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

    cplexroot="/opt/ibm/ILOG/CPLEX_Studio1261"

    cmake -Wdev --debug-trycompile \
        -D GA_DIR:STRING="$prefix/ga-5-4" \
        -D USE_PROGRESS_RANKS:BOOL=OFF \
        -D BOOST_ROOT:STRING="$prefix" \
        -D PETSC_DIR:STRING="/net/flophouse/files0/perksoft/petsc-3.7.2" \
        -D PETSC_ARCH:STRING='linux-gnu48-real-opt' \
        -D MPI_CXX_COMPILER:STRING="$prefix/bin/mpicxx" \
        -D MPI_C_COMPILER:STRING="$prefix/bin/mpicc" \
        -D MPIEXEC:STRING="$prefix/bin/mpiexec" \
        -D USE_CPLEX:BOOL=ON \
        -D CPLEX_ROOT_DIR:PATH="$cplexroot" \
        -D MPIEXEC_MAX_NUMPROCS:STRING="4" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=10 \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack" \
        $common_flags ..

elif [ $host == "flophouse44" ]; then

    # RHEL 5 with stock GNU 4.4 compilers

    prefix="/net/flophouse/files0/perksoft/linux64/openmpi44"
    PATH="${prefix}/bin:${PATH}"
    export PATH

    CC="gcc44"
    export CC
    CXX="g++44"
    export CXX
    CFLAGS="-pthread -Wall"
    export CFLAGS
    CXXFLAGS="-pthread -Wall"
    export CXXFLAGS

    cmake -Wdev --debug-trycompile \
        -D GA_DIR:STRING="$prefix" \
        -D USE_PROGRESS_RANKS:BOOL=OFF \
        -D BOOST_ROOT:STRING="$prefix" \
        -D PETSC_DIR:STRING="/net/flophouse/files0/perksoft/petsc-3.4.3" \
        -D PETSC_ARCH:STRING='arch-linux2-g++44-opt' \
        -D MPI_CXX_COMPILER:STRING="$prefix/bin/mpicxx" \
        -D MPI_C_COMPILER:STRING="$prefix/bin/mpicc" \
        -D MPIEXEC:STRING="$prefix/bin/mpiexec" \
        -D MPIEXEC_MAX_NUMPROCS:STRING="4" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=10 \
        -D USE_GLPK:BOOL=ON \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack" \
        $common_flags ..
    
elif [ $host == "pe10900" ]; then

    # Mac using GNU 4.8 and OpenMPI -- avoid using the system
    # compilers and MPI wrappers -- using MacPorts

    CC=/opt/local/bin/gcc
    export CC
    CXX=/opt/local/bin/g++
    export CXX

    prefix="/net/flophouse/files0/perksoft/macosx"
    cplexroot="/opt/ibm/ILOG/CPLEX_Studio1261/"

    cmake $options \
        -D GA_DIR:STRING="$prefix" \
        -D GA_EXTRA_LIBS:STRING="-lblas" \
        -D BOOST_ROOT:STRING='/opt/local' \
        -D PETSC_DIR:STRING="$prefix/../petsc-3.7.2" \
        -D PETSC_ARCH:STRING='arch-macosx-complex-opt' \
        -D MPI_CXX_COMPILER:STRING='/opt/local/bin/mpicxx' \
        -D MPI_C_COMPILER:STRING='/opt/local/bin/mpicc' \
        -D MPIEXEC:STRING='/opt/local/bin/mpiexec' \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=10 \
        -D USE_CPLEX:BOOL=ON \
        -D CPLEX_ROOT_DIR:PATH="$cplexroot" \
        -D USE_GLPK:BOOL=OFF \
        -D GLPK_ROOT_DIR:PATH="/opt/local" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack" \
        $common_flags ..

elif [ $host == "pe10900intel" ]; then

    # CMake really, really likes to use the wrong compiler. This
    # system uses MacPorts to supply a GNU 4.8 compiler. In order to
    # get GridPACK to build with the Intel compilers, the MacPorts
    # compilers need to be avoided. Do this as root:
    # 
    # port select gcc none

    prefix="/opt/intel/openmpi"
    PATH="$prefix/bin:$PATH" 
    RPATH="$prefix/lib:/opt/intel/lib"
    DYLD_LIBRARY_PATH="$RPATH"

    CC=icc
    CXX=icpc
    CFLAGS="-static-intel"
    CXXFLAGS="-static-intel"

    export PATH CC CXX CFLAGS CXXFLAGS RPATH DYLD_LIBRARY_PATH

    cmake -Wdev --debug-trycompile \
        -D GA_DIR:STRING="$prefix" \
        -D BOOST_ROOT:STRING="$prefix" \
        -D PETSC_DIR:STRING="$prefix/petsc-3.6.0" \
        -D PETSC_ARCH:STRING="arch-macosx-complex-opt" \
        -D MPI_CXX_COMPILER:STRING="$prefix/bin/mpicxx" \
        -D MPI_C_COMPILER:STRING="$prefix/bin/mpicc" \
        -D MPIEXEC:STRING="$prefix/bin/mpiexec" \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=10 \
        -D USE_GLPK:BOOL=OFF \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack" \
        -D CMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        ..


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

elif [ $host == "gridpackvm" ]; then

    prefix="$HOME/gridpack"
    root="/opt/ibm/ILOG/CPLEX_Studio1261/"
    cmake -Wno-dev --debug-try-compile \
	-D PETSC_DIR:STRING="$prefix/petsc-3.6.2" \
	-D PETSC_ARCH:STRING="arch-linux-real-opt" \
	-D GA_DIR:STRING="$prefix" \
	-D MPI_CXX_COMPILER:STRING="mpicxx" \
	-D MPI_C_COMPILER:STRING="mpicc" \
	-D MPIEXEC:STRING="mpiexec" \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=20 \
        -D USE_GLPK:BOOL=OFF \
        -D GLPK_ROOT_DIR:PATH="/opt/local" \
        -D USE_CPLEX:BOOL=ON \
        -D CPLEX_ROOT_DIR:PATH="$root" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack" \
	$common_flags ..

else

    echo "Unknown host: $host"
    exit 2
    
fi

