# -------------------------------------------------------------
# file: BuildGAExternalProject.cmake
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created October 12, 2018 by William A. Perkins
# Last Change: 2019-08-07 14:51:01 d3g096
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
    list(APPEND GA_OPTS
      -D BUILD_SHARED_LIBS:BOOL=YES
      -D CMAKE_POSITION_INDEPENDENT_CODE:BOOL=YES
      )
  else()
    list(APPEND GA_OPTS -D BUILD_SHARED_LIBS:BOOL=NO )
  endif()

  if (USE_PROGRES_RANKS)
    list(APPEND GA_OPTS -D GA_RUNTIME:STRING=MPI_PROGRESS_RANK)
  else ()
    list(APPEND GA_OPTS -D GA_RUNTIME:STRING=MPI_2SIDED)
  endif()

  include(ExternalProject)
  ExternalProject_Add(external_global_arrays
    PREFIX ${BUILD_DIR}/ga
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/ga
    INSTALL_DIR ${BIN_DIR}/ga
    CMAKE_ARGS -D ENABLE_BLAS:BOOL=NO
    -D ENABLE_FORTRAN:BOOL=NO
    -D ENABLE_CXX:BOOL=YES
    -D MPI_CXX_COMPILER:STRING=${MPI_CXX_COMPILER}
    -D MPI_C_COMPILER:STRING=${MPI_C_COMPILER}
    ${GA_OPTS}
    ${DEFAULT_CMAKE_FLAGS}
    -D CMAKE_INSTALL_PREFIX=${BIN_DIR}/ga
    BUILD_IN_SOURCE 0
    BUILD_COMMAND make 
    INSTALL_COMMAND make install  
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