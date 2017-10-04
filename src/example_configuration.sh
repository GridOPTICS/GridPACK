#! /bin/sh


# -------------------------------------------------------------
# handle command line options
# -------------------------------------------------------------
usage="$0 [-d|-r] [-s] [name]"

set -- `getopt rds $*`
if [ $? != 0 ]; then
    echo $usage >&2
    exit 2
fi

build="RelWithDebInfo"
shared="OFF"
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
        -s)
            shared="ON"
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
        -D BUILD_SHARED_LIBS:BOOL=$shared \
        -D CMAKE_BUILD_TYPE:STRING=$build \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
"

if [ $host == "flophouse" ]; then

    # RHEL 7 with GNU 4.8 compilers
    # The following are installed as packages:
    #   boost.x86_64                         1.53.0-26.el7
    #   boost-openmpi.x86_64                 1.53.0-26.el7
    #   (lots of other boost packages)
    #   

    prefix="/net/flophouse/files0/perksoft/linux64"
    PATH="${prefix}/bin:${PATH}"
    export PATH

    CC="/usr/bin/gcc"
    export CC
    CXX="/usr/bin/g++"
    export CXX

    cplexroot="/opt/ibm/ILOG/CPLEX_Studio1261"

    cmake -Wdev --debug-trycompile \
        -D GA_DIR:STRING="$prefix/ga-c++" \
        -D USE_PROGRESS_RANKS:BOOL=OFF \
        -D BOOST_ROOT:STRING="/usr" \
        -D PETSC_DIR:STRING="/net/flophouse/files0/perksoft/petsc-3.7.6" \
        -D PETSC_ARCH:STRING='linux-gnu48-real-opt' \
        -D MPI_CXX_COMPILER:STRING="/usr/lib64/openmpi/bin/mpicxx" \
        -D MPI_C_COMPILER:STRING="/usr/lib64/openmpi/bin/mpicc" \
        -D MPIEXEC:STRING="/usr/lib64/openmpi/bin/mpiexec" \
        -D USE_CPLEX:BOOL=OFF \
        -D CPLEX_ROOT_DIR:PATH="$cplexroot" \
        -D USE_GLPK:BOOL=ON \
        -D MPIEXEC_MAX_NUMPROCS:STRING="4" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=10 \
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
        -D USE_GLPK:BOOL=ON \
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

elif [ $host == "WE32673" ]; then

    # Mac using CLang 3.8 compilers and OpenMPI via MacPorts
    # The following MacPorts packages are installed:
    #   clang-3.8 @3.8.1_8+analyzer
    #   openmpi-clang38 @1.10.3_0+gcc6
    #   parmetis @4.0.3_3+openmpi
    #   boost @1.59.0_2+clang38+no_single+openmpi+python27
    #   global-arrays @5.5_10746+clang38+cxx+openmpi

    CC=/opt/local/bin/clang-mp-3.8
    export CC
    CXX=/opt/local/bin/clang++-mp-3.8
    export CXX

    prefix="/Users/d3g096/Projects/GridPACK"

    cmake $options \
        -D GA_DIR:STRING="/opt/local" \
        -D BOOST_ROOT:STRING="/opt/local" \
        -D PETSC_DIR:STRING="$prefix/petsc-3.7.5" \
        -D PETSC_ARCH:STRING="arch-macosx-clang-real-opt" \
        -D MPI_CXX_COMPILER:STRING='/opt/local/bin/mpicxx' \
        -D MPI_C_COMPILER:STRING='/opt/local/bin/mpicc' \
        -D MPIEXEC:STRING='/opt/local/bin/mpiexec' \
        -D MPIEXEC_MAX_NUMPROCS:STRING="4" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=10 \
        -D USE_CPLEX:BOOL=OFF \
        -D USE_GLPK:BOOL=ON \
        -D GLPK_ROOT_DIR:PATH="/opt/local" \
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

elif [ $host == "gridpackvm" ]; then

    prefix="$HOME/gridpack"

    CC=gcc
    CXX=g++
    export CC CXX

    cmake -Wno-dev --debug-try-compile \
	-D PETSC_DIR:STRING="$prefix/petsc-3.7.5" \
	-D PETSC_ARCH:STRING="arch-ubuntu-real-opt" \
	-D GA_DIR:STRING="$prefix" \
	-D MPI_CXX_COMPILER:STRING="mpicxx" \
	-D MPI_C_COMPILER:STRING="mpicc" \
	-D MPIEXEC:STRING="mpiexec" \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=20 \
        -D USE_GLPK:BOOL=ON \
        -D GLPK_ROOT_DIR:PATH="/usr" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix" \
	$common_flags ..

else

    echo "Unknown host: $host"
    exit 2
    
fi

