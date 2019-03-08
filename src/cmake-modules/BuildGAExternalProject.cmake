# -------------------------------------------------------------
# file: BuildGAExternalProject.cmake
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created October 12, 2018 by William A. Perkins
# Last Change: 2019-03-08 12:57:01 d3g096
# -------------------------------------------------------------

# -------------------------------------------------------------
# BuildGAExternalProject 
#
# This downloads and builds a recent version of the Global Arrays
# library.  The following variables are set to mimic find_package(GA):
#
#  GA_FOUND        - system has GA
#  GA_DIR          - root directory of GA installation
#  GA_INCLUDE_DIRS - include directories for GA
#  GA_LIBRARIES    - libraries for GA
# -------------------------------------------------------------
function(BuildGAExternalProject)
  # FIXME: allow options for GA communication options: ts, pr, openib?
  # and allow for extra libraries.


  set(GA_OPTS "")
  set(GA_LIBS "")

  if(BUILD_SHARED_LIBS)
    list(APPEND GA_OPTS --enable-shared=yes --enable-static=no)
  else()
    list(APPEND GA_OPTS --enable-shared=no --enable-static=yes)
  endif()

  if(GA_INFINIBAND) 
    list(APPEND GA_OPTS "--with-openib")
    set(GA_LIBS  "LIBS=-libverbs")
  else()
    list(APPEND GA_OPTS "--with-mpi-ts")
  endif()

  include(ExternalProject)
  ExternalProject_Add(external_global_arrays
    URL https://github.com/GlobalArrays/ga/releases/download/v5.7/ga-5.7.tar.gz
    URL_HASH MD5=bb9a441a6b4fbb8b52b58c2d3f4cd07f
    DOWNLOAD_DIR ${DOWNLOAD_DIR}
    CONFIGURE_COMMAND sh ./configure 
    --prefix=${BIN_DIR}/ga 
    --enable-cxx 
    --disable-f77 
    --enable-i4 
    --with-mpi 
    --without-blas 
    --without-lapack 
    --without-scalapack 
    ${GA_OPTS}
    CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    MPICC=${MPI_C_COMPILER} MPICXX=${MPI_CXX_COMPILER}
    ${GA_LIBS}
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make 
    INSTALL_COMMAND make install  
    PREFIX ${BUILD_DIR}/ga
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${BIN_DIR}/ga ${DEFAULT_CMAKE_FLAGS}
    INSTALL_DIR ${BIN_DIR}/ga
    )

  ExternalProject_Get_Property(external_global_arrays INSTALL_DIR)

  set(GA_DIR ${BIN_DIR}/ga PARENT_SCOPE)
  set(GA_DIR "${INSTALL_DIR}" PARENT_SCOPE)
  message(STATUS "${INSTALL_DIR}")
  message(STATUS "GA_DIR=${GA_DIR}")

  set(GA_DIR "${INSTALL_DIR}" PARENT_SCOPE)

  set(GA_FOUND TRUE PARENT_SCOPE)

  set(GA_INCLUDE_DIRS ${INSTALL_DIR}/include PARENT_SCOPE)

  set(GA_LIBRARIES
    ${INSTALL_DIR}/lib/libga++${LIB_SUFFIX}
    ${INSTALL_DIR}/lib/libga${LIB_SUFFIX}
    ${INSTALL_DIR}/lib/libarmci${LIB_SUFFIX}
    ${GA_LIBS} PARENT_SCOPE)
  
  # make sure there is a slash at the end of the path
  install(
    DIRECTORY ${BIN_DIR}/ga/include/
    DESTINATION include
    FILES_MATCHING PATTERN "*"
    )
  install(
    DIRECTORY ${BIN_DIR}/ga/lib/
    DESTINATION lib
    FILES_MATCHING PATTERN "*"
    )
  install(
    DIRECTORY ${BIN_DIR}/ga/bin/
    DESTINATION bin
    FILES_MATCHING PATTERN "*"
    )

endfunction(BuildGAExternalProject)