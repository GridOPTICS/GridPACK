#! /bin/sh


# -------------------------------------------------------------
# handle command line options
# -------------------------------------------------------------
usage="$0 [-d|-r] [-s] [-G] [-B] [name]"

options=`getopt rdsGB $*`
if [ $? != 0 ]; then
    echo $usage >&2
    exit 2
fi
set -- $options

build="RelWithDebInfo"
shared="FALSE"
buildGA="FALSE"
buildBoost="FALSE"
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

elif [ $host == "we32673" ]; then

    # Mac using CLang 14.0 compilers and MPICH via MacPorts
    # The following MacPorts packages are installed:
    #   clang-14 @14.0.6_0+analyzer+libstdcxx
    #   mpich-clang14 @4.0.2_0+gcc12
    #   boost176 @1.76.0_3+clang14+mpich+no_single+no_static+python310
    # Global Arrays 5.7 built by hand
    # PETSc 3.8.4 w/ ParMETIS, SuperLU, etc., built by hand
    # Need to make sure the compiler set and MPI are selected, i.e.
    #   sudo port select clang mp-clang-6.0
    #   sudo port select mpi mpich-clang60-fortran
    # Cannot use PETSc < 3.8.0

    CC=/opt/local/bin/clang
    export CC
    CXX=/opt/local/bin/clang++
    export CC CXX 

    prefix="/Users/d3g096/Projects/GridPACK"
    pdir="$prefix/petsc.gitlab"

    if [ "$shared"x = "ON"x ]; then
        parch="darwin-mpich-clang-real-opt-3.16"
        parch="darwin-mpich-clang-complex-opt-3.16"
    else
        parch="darwin-mpich-clang-complex-opt-3.16-static"
    fi
    
    cmake $options \
        -D GA_DIR:STRING="$prefix/gridpack-install" \
        -D BOOST_ROOT:STRING="/opt/local/libexec/boost/1.76" \
        -D Boost_NO_BOOST_CMAKE:BOOL=TRUE \
        -D PETSC_DIR:PATH="$pdir" \
        -D PETSC_ARCH:STRING="$parch" \
        -D MPI_CXX_COMPILER:STRING='/opt/local/bin/mpicxx' \
        -D MPI_C_COMPILER:STRING='/opt/local/bin/mpicc' \
        -D MPIEXEC:STRING='/opt/local/bin/mpiexec' \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=60 \
        -D USE_CPLEX:BOOL=OFF \
        -D USE_GLPK:BOOL=ON \
        -D GLPK_ROOT_DIR:PATH="/opt/local" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack-install" \
        $common_flags ..

elif [ $host == "WE30729" ]; then

    # Macbook using CLang 6.0 compilers and MPICH via MacPorts
    # Pretty much the same as WE32673

    CC=/opt/local/bin/clang
    export CC
    CXX=/opt/local/bin/clang++
    export CXX

    prefix="$HOME/Projects/GridPACK"

    if [ "$shared"x = "ON"x ]; then
        pdir="$prefix/petsc-3.10.5" 
        parch="macosx-complex-c-shared" 
    else
        pdir="$prefix/petsc-3.8.4"
        parch="arch-macosx-clang-real-opt"
    fi
    cmake $options \
        -D BOOST_ROOT:STRING="/opt/local" \
        -D PETSC_DIR:PATH="$pdir" \
        -D PETSC_ARCH:STRING="$parch" \
        -D MPI_CXX_COMPILER:STRING='/opt/local/bin/mpicxx' \
        -D MPI_C_COMPILER:STRING='/opt/local/bin/mpicc' \
        -D MPIEXEC:STRING='/opt/local/bin/mpiexec' \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=60 \
        -D USE_CPLEX:BOOL=OFF \
        -D USE_GLPK:BOOL=ON \
        -D GLPK_ROOT_DIR:PATH="/opt/local" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack-hadrec" \
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
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc.gitlab" \
    #      -D PETSC_ARCH:STRING="ubuntu-complex-shared-3.16.6" \
    #      -D USE_OLD_PETSC:BOOL=OFF \

    # Custom built 3.16, real:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc.gitlab" \
    #      -D PETSC_ARCH:STRING="ubuntu-real-shared-3.16.6" \
    #      -D USE_OLD_PETSC:BOOL=OFF \

    # System PETSc package:
    #      -D PETSC_DIR:STRING="/usr/lib/petsc" \
    #      -D PETSC_ARCH:STRING="" \
    #      -D USE_OLD_PETSC:BOOL=OFF \

    # Custom built 3.12.5, real:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.12.5" \
    #      -D PETSC_ARCH:STRING="ubuntu-real-shared-3.12" \
    #      -D USE_OLD_PETSC:BOOL=OFF \
    
    # Custom built 3.10.5, real:
    #      -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc-3.10.5" \
    #      -D PETSC_ARCH:STRING="ubuntu-real-shared-3.10" \
    #      -D USE_OLD_PETSC:BOOL=OFF \
    
    
    prefix="$HOME/Projects/GridPakLDRD/gridpack-install"
    cmake -Wdev --debug-trycompile \
          -D PETSC_DIR:STRING="/home/d3g096/Projects/GridPakLDRD/petsc.gitlab" \
          -D PETSC_ARCH:STRING="ubuntu-real-shared-3.16.6" \
          -D USE_OLD_PETSC:BOOL=OFF \
          -D BOOST_ROOT:PATH="/usr" \
          -D Boost_NO_BOOST_CMAKE:BOOL=TRUE \
          -D PARMETIS_DIR:PATH="/usr" \
          -D MPI_CXX_COMPILER:STRING="mpicxx" \
          -D MPI_C_COMPILER:STRING="mpicc" \
          -D MPIEXEC:STRING="mpiexec" \
          -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
          -D GRIDPACK_TEST_TIMEOUT:STRING=60 \
          -D USE_GLPK:BOOL=OFF \
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

