rm -rf CMake*
CC=gcc
CXX=g++
CXXFLAGS=-lrt
export CC CXX CXXFLAGS

cmake \
    -D PETSC_DIR:STRING="/home/huan495/petsc-3.7.7" \
    -D PARMETIS_DIR:PATH="/usr" \
    -D GA_EXTRA_LIBS:STRING="-lscalapack-openmpi -lblacsCinit-openmpi -lblacs-openmpi -llapack -lblas -lgfortran -lrt" \
    -D MPI_CXX_COMPILER:STRING="mpicxx" \
    -D MPI_C_COMPILER:STRING="mpicc" \
    -D MPIEXEC:STRING="mpiexec" \
    -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
    -D GRIDPACK_TEST_TIMEOUT:STRING=60 \
    -D USE_GLPK:BOOL=ON \
    -D GLPK_ROOT_DIR:PATH="/usr" \
    -D BUILD_SHARED_LIBS:BOOL=OFF \
    -D CMAKE_INSTALL_PREFIX:PATH="$HOME/gridpack" \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
    -D HELICS_INSTALL_DIR:PATH="/home/huan495/HELICS/install-new" \
    -D ZEROMQ_INSTALL_DIR:PATH="/home/huan495/zmq-install" \
    CXXFLAGS=-lrt \
    ..
