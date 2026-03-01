#! /bin/sh


# -------------------------------------------------------------
# handle command line options
# -------------------------------------------------------------
usage="$0 [-d|-r] [-s] [-G] [-B] [-c] [name]"

options=`getopt crdsGB $*`
if [ $? != 0 ]; then
    echo $usage >&2
    exit 2
fi
set -- $options

build="RelWithDebInfo"
shared="FALSE"
buildGA="FALSE"
buildBoost="FALSE"
initcomm="OFF"
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
        -G)
            buildGA="ON"
            shift
            ;;
        -B)
            buildBoost="ON"
            shift
            ;;
        -c)
            initcomm="ON"
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
        -D BUILD_GA:BOOL=$buildGA \
        -D BUILD_BOOST:BOOL=$buildBoost \
        -D BUILD_SHARED_LIBS:BOOL=$shared \
        -D CMAKE_BUILD_TYPE:STRING=$build \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        -D ENABLE_ENVIRONMENT_FROM_COMM:BOOL=$initcomm \
"

if [ $host == "briareus" ]; then

    # Using GNU 4.9 and OpenMPI modules

    prefix="/files0/gridpack.gnu"
    PATH="${prefix}/bin:${PATH}"
    export PATH

    CC="/share/apps/gcc/4.9.2/bin/gcc"
    export CC
    CXX="/share/apps/gcc/4.9.2/bin/g++"
    export CXX

    cmake -Wdev --debug-trycompile \
        -D USE_PROGRESS_RANKS:BOOL=OFF \
        -D GA_DIR:PATH="$prefix" \
        -D BOOST_ROOT:STRING="$prefix" \
        -D PETSC_DIR:STRING="/files0/petsc-3.7.5" \
        -D PETSC_ARCH:STRING="gridpack-gnu-openmpi-real" \
        -D MPI_CXX_COMPILER:STRING="mpicxx" \
        -D MPI_C_COMPILER:STRING="mpicc" \
        -D MPIEXEC:STRING="mpiexec" \
        -D USE_CPLEX:BOOL=OFF \
        -D USE_GLPK:BOOL=OFF \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=10 \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix" \
        $common_flags ..

elif [ $host == "WE39945" ]; then

    # Macbook Pro using CLang 16.0 compilers and MPICH via MacPorts
    # The following MacPorts packages are installed:
    #   clang-16 @16.0.6_6+analyzer
    #   mpich-clang16 @4.2.3_0+gcc14
    # Global Arrays 5.7 built by hand
    # PETSc 3.8.4 w/ ParMETIS, SuperLU, etc., built by hand
    # Boost 1.81 w/ MPI build by hand
    # Need to make sure the compiler set and MPI are selected, i.e.
    #   sudo port select --set mpi mpich-clang16-fortran
    #   sudo port select clang mp-clang-6.0

    CC=clang
    export CC
    CXX=clang++
    export CC CXX 

    prefix="$HOME/Projects/GridPACK/install"

    cmake $options \
        --graphviz=GridPACK.dot \
        -D GA_DIR:STRING="$prefix" \
        -D Boost_ROOT:STRING="/opt/local/libexec/boost/1.81" \
        -D PETSC_DIR:PATH="/Users/d3g096/Projects/GridPACK/src/petsc-3.23.6" \
        -D PETSC_ARCH:STRING="mpich-clang-real-shared" \
        -D MPI_CXX_COMPILER:STRING='/opt/local/bin/mpicxx' \
        -D MPI_C_COMPILER:STRING='/opt/local/bin/mpicc' \
        -D MPIEXEC:STRING='/opt/local/bin/mpiexec' \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=120 \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix" \
        $common_flags ..

elif [ $host == "constance" ]; then

    CC=`which gcc`
    CXX=`which g++`
    CFLAGS="-pthread"
    CXXFLAGS="-pthread"

    export CC CXX CFLAGS CXXFLAGS

    prefix="/people/d3g096/gridpack"
    cmake $options \
        -D GA_DIR:STRING="$prefix/ga-5-4" \
	-D GA_EXTRA_LIBS:STRING="-libverbs" \
	-D BOOST_ROOT:STRING="$prefix" \
	-D PETSC_DIR:STRING="$prefix/petsc-3.7.3" \
	-D PETSC_ARCH:STRING='constance-gnu48-complex-opt' \
	-D MPI_CXX_COMPILER:STRING=`which mpicxx` \
	-D MPI_C_COMPILER:STRING=`which mpicc` \
	-D MPIEXEC:STRING=`which mpiexec` \
        -D GRIDPACK_TEST_TIMEOUT:STRING=10 \
	$common_flags ..

elif [ $host == "constance-hadrec" ]; then

    # always build shared
    # always build GA

    CC=`which gcc`
    CXX=`which g++`
    CFLAGS="-pthread"
    CXXFLAGS="-pthread"
    export CC CXX CFLAGS CXXFLAGS

    prefix="/pic/projects/gripdack/hadrec/gridpack"
    cmake $options \
	-D BOOST_ROOT:STRING="$prefix" \
	-D PETSC_DIR:STRING="$prefix/petsc-3.9.4" \
	-D PETSC_ARCH:STRING='linux-openmpi-gnu-cxx-complex-opt' \
	-D MPI_CXX_COMPILER:STRING=`which mpicxx` \
	-D MPI_C_COMPILER:STRING=`which mpicc` \
	-D MPIEXEC:STRING=`which mpiexec` \
        -D GRIDPACK_TEST_TIMEOUT:STRING=10 \
        -D BUILD_GA:BOOL=ON \
        -D BUILD_SHARED_LIBS:BOOL=ON \
        -D CMAKE_BUILD_TYPE:STRING=$build \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=FALSE \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack" \
	..

elif [ $host == "tlaloc" ]; then

    # Ubuntu 20 with as many system packages as possible: GNU
    # compilers 9.4.0, OpenMPI 4.0.3, Boost 1.71.0, PETSc 3.12,
    # ParMETIS 4.0.3.  
    
    CC=gcc
    CXX=g++
    CFLAGS=-pthread
    CXXFLAGS=-pthread
    export CC CXX CFLAGS CXXFLAGS

    # Custom built 3.16, complex:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.16.6" \
    #      -D PETSC_ARCH:STRING="ubuntu-complex-shared" \

    # Custom built 3.16, real:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.16.6" \
    #      -D PETSC_ARCH:STRING="ubuntu-real-shared" \

    # Custom built 3.16, complex, w/o superlu_dist:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.16.6" \
    #      -D PETSC_ARCH:STRING="ubuntu-complex-shared-mumps" \

    # Custom built 3.19, complex:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.19.4" \
    #      -D PETSC_ARCH:STRING="ubuntu-complex-shared" \

    # Custom built 3.19, real:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.19.4" \
    #      -D PETSC_ARCH:STRING="ubuntu-real-shared" \

    # Custom built 3.19, complex, w/o superlu_dist:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.19.4" \
    #      -D PETSC_ARCH:STRING="ubuntu-complex-shared-mumps" \

    # Custom built 3.20, complex:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.20.6" \
    #      -D PETSC_ARCH:STRING="ubuntu-complex-shared" \

    # Custom built 3.20, real:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.20.6" \
    #      -D PETSC_ARCH:STRING="ubuntu-real-shared" \

    # Custom built 3.21, complex:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.21.6" \
    #      -D PETSC_ARCH:STRING="ubuntu-complex-shared-debug" \

    # Custom built 3.21, complex:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.21.6" \
    #      -D PETSC_ARCH:STRING="ubuntu-real-shared-debug" \

    if [ -z "$GRIDPACK_DIR" ]; then
        prefix="$HOME/Projects/GridPACK-Wind/gridpack-install"
    else
        prefix="$GRIDPACK_DIR"
    fi
    
    cmake -Wdev --debug-trycompile \
          --graphviz=GridPACK.dot \
          -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.20.6" \
          -D PETSC_ARCH:STRING="ubuntu-complex-shared" \
          -D Boost_ROOT:PATH="/usr" \
          -D PARMETIS_DIR:PATH="/usr" \
          -D MPI_CXX_COMPILER:STRING="mpicxx" \
          -D MPI_C_COMPILER:STRING="mpicc" \
          -D MPIEXEC:STRING="mpiexec" \
          -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
          -D GRIDPACK_TEST_TIMEOUT:STRING=60 \
          -D USE_GLPK:BOOL=OFF \
          -D GA_DIR="$HOME/Projects/ExaLearn/ga-install" \
          -D BUILD_SHARED_LIBS:BOOL=ON \
          -D CMAKE_INSTALL_PREFIX:PATH="$prefix" \
          -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
          $common_flags ..

elif [ $host == "gridpackvm" ]; then

    prefix="$HOME/gridpack"

    CC=gcc
    CXX=g++
    export CC CXX

    cmake $options \
	-D PETSC_DIR:STRING="/usr/lib/petsc" \
	-D PARMETIS_DIR:PATH="/usr" \
	-D GA_EXTRA_LIBS:STRING="-lscalapack-openmpi -lblacsCinit-openmpi -lblacs-openmpi -llapack -lblas -lgfortran" \
	-D MPI_CXX_COMPILER:STRING="mpicxx" \
	-D MPI_C_COMPILER:STRING="mpicc" \
	-D MPIEXEC:STRING="mpiexec" \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=20 \
        -D USE_GLPK:BOOL=ON \
        -D GLPK_ROOT_DIR:PATH="/usr" \
        -D CMAKE_INSTALL_PREFIX:PATH="$HOME/gridpack/gridpack_install" \
	-D HELICS_INSTALL_DIR:PATH="/home/huan495/HELICS/install-test" \
	$common_flags ..

elif [ $host == "debianvm" ]; then

    prefix="$HOME/gridpack"

    CC=gcc
    CXX=g++
    CFLAGS=-pthread
    CXXFLAGS=-pthread
    export CC CXX CFLAGS CXXFLAGS

    cmake $options \
	-D PETSC_DIR:STRING="/usr/lib/petsc" \
	-D PARMETIS_DIR:PATH="/usr" \
	-D GA_EXTRA_LIBS:STRING="-lscalapack-openmpi -lblacs-openmpi -llapack -lblas -lgfortran" \
	-D MPI_CXX_COMPILER:STRING="mpicxx" \
	-D MPI_C_COMPILER:STRING="mpicc" \
	-D MPIEXEC:STRING="mpiexec" \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=30 \
        -D USE_GLPK:BOOL=ON \
        -D GLPK_ROOT_DIR:PATH="/usr" \
        -D BUILD_SHARED_LIBS:BOOL=OFF \
        -D CMAKE_INSTALL_PREFIX:PATH="/usr" \
	$common_flags ..

else

    echo "Unknown host: $host"
    exit 2
    
fi

