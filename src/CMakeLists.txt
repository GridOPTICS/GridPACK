#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be
#     found
#     in the LICENSE file in the top level directory of this distribution.
#
# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May  3, 2013 by William A. Perkins
# Last Change: 2019-11-14 06:49:07 d3g096
# -------------------------------------------------------------

# 
cmake_minimum_required (VERSION 3.3.0)

# This allows the use of LOCATION(target)
if (POLICY CMP0026)
  cmake_policy(SET CMP0026 OLD)
endif (POLICY CMP0026)


project (GridPACK)

set (GridPACK_VERSION_MAJOR 3)
set (GridPACK_VERSION_MINOR 3)
set (GridPACK_VERSION_PATCH 0)

# This is a C++ project exclusively
enable_language(CXX)

# where to look for special .cmake files
list (APPEND CMAKE_MODULE_PATH 
    "${GridPACK_SOURCE_DIR}/cmake-modules"
    "${GridPACK_SOURCE_DIR}/../cmake-jedbrown"
)
include(GridPACK)

# add user configuration options
option (USE_PROGRESS_RANKS "Enable progress ranks in GridPACK" OFF)
if (USE_PROGRESS_RANKS)
  add_definitions (-DUSE_PROGRESS_RANKS=1)
endif()

# add GOSS directory
option (GOSS_DIR "Point to directory with GOSS files" OFF)
if (GOSS_DIR)
  add_definitions (-DUSE_GOSS=1) 
  include_directories(AFTER ${GOSS_DIR}/include/activemq-cpp-3.8.4)
  set(GOSS_INCLUDE "${GOSS_DIR}/include/activemq-cpp-3.8.4")
  set (GOSS_LIBRARY "${GOSS_DIR}/lib/libactivemq-cpp.a")
  message(STATUS "Set GOSS library at ${GOSS_LIBRARY}")
endif()
if (APR_DIR)
  include_directories(AFTER ${APR_DIR}/include/apr-1)
  set (APR_INCLUDE "${APR_DIR}/include/apr-1")
  set (APR_LIBRARY "${APR_DIR}/lib/libapr-1.a")
  message(STATUS "Set APR library at ${APR_LIBRARY}")
endif()

# check to see if rt library exists. If so, add -lrt to GA_EXTRA_LIBS
find_library(LIBRT rt PATHS /lib /usr/lib /usr/local/lib)
if (LIBRT)
  if (GA_EXTRA_LIBS)
    set (GA_EXTRA_LIBS "${GA_EXTRA_LIBS} ${LIBRT}")
  else()
    set (GA_EXTRA_LIBS ${LIBRT})
  endif()
endif()

# add HELICS directory
option (HELICS_INSTALL_DIR "Point to directory with installed HELICS files" OFF)
if (HELICS_INSTALL_DIR AND ZEROMQ_INSTALL_DIR)
  message(STATUS "Building GridPACK with HELICS library")
  include_directories(AFTER ${HELICS_INSTALL_DIR}/include)
  if (EXISTS ${HELICS_INSTALL_DIR}/lib)
    message(STATUS "Found HELICS library at ${HELICS_INSTALL_DIR}lib")
    set (HELICS_LIBRARY "${HELICS_INSTALL_DIR}/lib/libhelics-static.a")
    set (JSON_LIBRARY "${HELICS_INSTALL_DIR}/lib/libjsoncpp.a")
  elseif (EXISTS ${HELICS_INSTALL_DIR}/lib64)
    message(STATUS "Found HELICS library at ${HELICS_INSTALL_DIR}lib64")
    set (HELICS_LIBRARY "${HELICS_INSTALL_DIR}/lib64/libhelics-static.a")
    set (JSON_LIBRARY "${HELICS_INSTALL_DIR}/lib64/libjsoncpp.a")
  endif()
  if (EXISTS ${ZEROMQ_INSTALL_DIR}/lib/libzmq.a)
    message(STATUS "Found ZeroMQ library at ${ZEROMQ_INSTALL_DIR}lib")
    set (ZEROMQ_LIBRARY "${ZEROMQ_INSTALL_DIR}/lib/libzmq.a")
    if (EXISTS ${ZEROMQ_INSTALL_DIR}/lib/libsodium.a)
      set (SODIUM_LIBRARY "${ZEROMQ_INSTALL_DIR}/lib/libsodium.a")
    endif()
  endif()
  if (HELICS_LIBRARY AND ZEROMQ_LIBRARY)
    add_definitions(-DUSE_HELICS=1)
  else()
    message(STATUS "HELICS is not configured. Not all libraries found")
  endif()
endif()

# suppress any checks to see if test codes run. Only check for compilation.
# use for cross-compilation situations
option (CHECK_COMPILATION_ONLY "Check compilation only" OFF)

# should cplex or gurobi be used
option(USE_CPLEX
  "Use optimization software CPLEX libraries." OFF)

option(USE_GUROBI
  "Use optimization software GUROBI libraries." OFF)

option(USE_GLPK
  "Use GNU linear programming kit optimization software library." OFF)

# BUILD_SHARED_LIBS is a CMake option and is understood by the CMake build
# system, independently of declaring it here
option(BUILD_SHARED_LIBS
  "Attempt to build all libraries as shared" OFF)


# turn testing on/off
option(GRIDPACK_ENABLE_TESTS "Enable build of unit tests" ON)
if (GRIDPACK_ENABLE_TESTS)
  enable_testing()
endif()

# -------------------------------------------------------------
# External project settings, generally
# -------------------------------------------------------------

# A set of common places to put the external projects (if used)

set(DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Downloads)
set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/build)
set(BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/libs)
set(GRIDPACK_TEST_LDPATH "")
if(BUILD_SHARED_LIBS)
  set(LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(LIB_PREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
  set(LIB_TYPE SHARED)
else()
  set(LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(LIB_PREFIX ${CMAKE_STATIC_LIBRARY_PREFIX})
  set(LIB_TYPE STATIC)
endif()

# make a bogus target so that all libraries depend on; external
# projects are added as dependencies to this target so external
# libraries are (hopefully) built before any GridPACK code

add_custom_target(external_build
  COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/external_build.stamp
  )

# -------------------------------------------------------------
# MPI compiler
# -------------------------------------------------------------
message(STATUS "Checking MPI ...")
find_package (MPI REQUIRED)
include_directories(AFTER ${MPI_CXX_INCLUDE_PATH})

if (ENABLE_CRAY_BUILD)
  set (CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

# This is here for cmake version 2.6 
if (NOT MPI_LIBRARY OR NOT MPI_EXTRA_LIBRARY)
# Punt if MPI_LIBRARY or MPI_EXTRA_LIBRARY not found
  set(MPI_CXX_LIBRARIES "")
else()
  if (NOT MPI_CXX_LIBRARIES)
   set(MPI_CXX_LIBRARIES 
     ${MPI_LIBRARY} 
     ${MPI_EXTRA_LIBRARY}
   )
  endif()
endif()
message(STATUS "MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")

# Perkins found out that this was exactly the wrong thing to do:
# 
# set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
# set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
# set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

# -------------------------------------------------------------
# Compiler specific options
# -------------------------------------------------------------

# FIXME: remove annoying warnings (these warnings probably identify
# real problems, but there are so many of them they obscure errors
# that are causing identifiable problems

if (CMAKE_COMPILER_IS_GNUCXX) 
  add_definitions(
    -pedantic
    -Wno-write-strings
    -Wno-long-long 
    -Wno-sign-compare 
    -Wno-unused-variable
    -Wno-unused-but-set-variable
    -Wno-maybe-uninitialized
    )
endif()

# using regular Clang or AppleClang
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_definitions(
    -Wno-c++11-extensions
    -Wno-c++11-compat-deprecated-writable-strings
    -Wno-format-security
    -DGRIDPACK_AVOID_APPLECLANG_MPI_PROBLEMS=1
    )
endif()


# When built as shared libraries, and installed in non-system dirs, deal with RPATH
if(BUILD_SHARED_LIBS)
  if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_MACOSX_RPATH 1)
  endif()
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES ${CMAKE_INSTALL_PREFIX}/lib isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
    set(CMAKE_SKIP_BUILD_RPATH  FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  else()
    set(CMAKE_SKIP_BUILD_RPATH  TRUE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
  endif()
endif()


# -------------------------------------------------------------
# Boost
#
# Find the Boost library. It needs to use the same compiler as
# MPI. How can I check that?
# -------------------------------------------------------------

option(BOOST_BUILD_PYTHON "Include Boost::Python in the Boost build" OFF)

# These are the Boost libraries required here
  
set(BOOST_BUILD_LIBS "")
list(APPEND BOOST_BUILD_LIBS mpi serialization random filesystem system)

if(BUILD_BOOST)

  if(BOOST_BUILD_PYTHON)
    list(APPEND BOOST_BUILD_LIBS python)
  endif()

  include(BuildBoostExternalProject)
  BuildBoostExternalProject("${BOOST_BUILD_LIBS}")
  add_dependencies(external_build external_boost)
  # FIXME: Build type (probably only important for Windoze)
  list(APPEND GRIDPACK_TEST_LDPATH ${Boost_LIBRARY_DIR_DEBUG} ${Boost_LIBRARY_DIR_RELEASE} )
else(BUILD_BOOST)
  message(STATUS "Checking Boost ...")
  find_package(Boost 1.53
    COMPONENTS ${BOOST_BUILD_LIBS}
    REQUIRED
    )
endif(BUILD_BOOST)
message(STATUS "Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES=${Boost_LIBRARIES}")

include_directories(AFTER ${Boost_INCLUDE_DIRS})

# -------------------------------------------------------------
# PETSc
#
# We need a very recent version of PETSc. It needs to parallel and
# compiled with C++ complex support.
# -------------------------------------------------------------
message(STATUS "Checking PETSc ...")
find_package (PETSc REQUIRED)
message(STATUS "Using PETSc version ${PETSC_VERSION}")

if (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/conf/PETScConfig.cmake")
  include("${PETSC_DIR}/${PETSC_ARCH}/conf/PETScConfig.cmake")
elseif (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/PETScConfig.cmake")
  include("${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/PETScConfig.cmake")
elseif (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/PETScBuildInternal.cmake")
  include("${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/PETScBuildInternal.cmake")
else()
  message(FATAL_ERROR "PETSc found, but CMake configuration for PETSc installation not found?")
endif()


# checks 

if (NOT PETSC_HAVE_MPI)
  message(FATAL_ERROR "PETSc installation is not parallel (--with-mpi=1)")
endif()

# PETSc can be compiled for double precsion (--with-precision=double)
# complex (--with-scalar-type=complex) or real
# (--with-scalar-type=real).  This is to determine what that
# underlying type is.

if (PETSC_USE_REAL_DOUBLE) 
  message(STATUS "PETSc installation is double precision (--with-precision=double) -- good")
else()
  message(FATAL_ERROR "PETSc installation is not double precision (--with-precision=double)")
endif()

if (PETSC_USE_COMPLEX) 
  message (STATUS "PETSc installation uses complex type (--with-scalar-type=complex)")
else()
  message (STATUS "PETSc installation uses real type (--with-scalar-type=real)")
endif()

# This is probably redundant
if (PETSC_CLANGUAGE_Cxx)
  message (STATUS "PETSc installation uses C++ (--with-clanguage=c++) -- we can work with that.")
endif()

if (PETSC_HAVE_SUPERLU_DIST) 
  set(GRIDPACK_MATSOLVER_PKG "superlu_dist")
  message(STATUS "PETSc parallel LU linear solver will be from SuperLU_dist")
elseif (PETSC_HAVE_MUMPS)
  set(GRIDPACK_MATSOLVER_PKG "mumps")
  message(STATUS "PETSc parallel LU linear solver will be from MUMPS")
else()
  set(GRIDPACK_MATSOLVER_PKG "petsc")
  message(WARNING "Installed PETSc does not have a parallel direct solver, expect unit test failures")
endif()

# -------------------------------------------------------------
# ParMETIS
# 
# Only versions > 4.0 have been used
# -------------------------------------------------------------
message(STATUS "Checking ParMETIS ...")

# If ParMETIS is included in PETSc, use that installation. This will
# override any specification on the command line.  
if (PETSC_HAVE_PARMETIS) 
  set(PARMETIS_DIR "${PETSC_DIR}/${PETSC_ARCH}")
endif()
find_package(ParMETIS REQUIRED)


# -------------------------------------------------------------
# GlobalArrays
# Build or Find
# -------------------------------------------------------------
option(BUILD_GA "Download and build Global Arrays" OFF)
option(GA_INFINIBAND "Try to configure GA for an Infiniband cluster" OFF)

if(BUILD_GA)
  include(BuildGAExternalProject)
  BuildGAExternalProject()
  add_dependencies(external_build external_global_arrays)
  list(APPEND GRIDPACK_TEST_LDPATH ${GA_DIR}/lib)
else(BUILD_GA)
  message(STATUS "Checking GA ...")
  find_package(GA REQUIRED)
endif(BUILD_GA)
message(STATUS "GA_INCLUDE_DIRS=${GA_INCLUDE_DIRS}")
message(STATUS "GA_LIBRARIES=${GA_LIBRARIES}")

# -------------------------------------------------------------
#  Add optimization libraries, if requested
# -------------------------------------------------------------
# add the CPLEX library?
if( USE_CPLEX)
  message(STATUS "Checking cplex ...")
  find_package(CPLEX REQUIRED)
endif( USE_CPLEX)

# add the GUROBI library?
IF( USE_GUROBI)
  message(STATUS "Checking Gurobi ...")
  find_package(GUROBI REQUIRED)
endif (USE_GUROBI)

if(USE_GLPK)
  message(STATUS "Checking glpk ...")
  find_package(GLPK REQUIRED)
endif(USE_GLPK)

#message(STATUS "Checking concert ...")
#find_package(CONCERT REQUIRED)

# -------------------------------------------------------------
# Doxygen
# -------------------------------------------------------------
find_package(Doxygen)
if(DOXYGEN_FOUND)

    # The custom target "doc" generates "user" documentation. This is
    # essentially the public interface of library classes.
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
      )

    

    # The custom target "devdoc" generates "developer"
    # documentation. This has all the nitty gritty details and code
    # that only people who work on GridPACK care about.
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.dev.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.dev @ONLY)
    add_custom_target(devdoc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.dev
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating devoloper documentation with Doxygen" VERBATIM
      )
endif(DOXYGEN_FOUND)

# -------------------------------------------------------------
# Options for all Programs/Libraries
# -------------------------------------------------------------

set(GRIDPACK_LIB_LINK_LIBRARIES OFF
  CACHE BOOL "Do we need library link libraries?")

if (BUILD_SHARED_LIBS)
  if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(GRIDPACK_LIB_LINK_LIBRARIES ON)
  endif()
endif()

set(GRIDPACK_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(GRIDPACK_DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/applications/data_sets)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (ENABLE_CRAY_BUILD)
  if (NOT ${MPI_LIBRARY} OR NOT ${MPI_EXTRA_LIBRARY})
    set(MPI_LIBRARY "")
    set(MPI_EXTRA_LIBRARY "")
  endif()
endif()

# If external packages are built as shared libraries, the library
# search path needs to be set so unit tests can find them when run.
# Try setting the LIBRARY_PATH environment variable to include the
# directories from external packages.

if(GRIDPACK_TEST_LDPATH)
  string(REPLACE ";" ":" tmp "${GRIDPACK_TEST_LDPATH}")
  set(GRIDPACK_TEST_LDPATH "${tmp}")
  if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    if(DEFINED ENV{DYLD_LIBRARY_PATH})
      set(GRIDPACK_TEST_LDPATH "${GRIDPACK_TEST_LDPATH}:$ENV{DYLD_LIBRARY_PATH}")
    endif()
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    if(DEFINED ENV{LD_LIBRARY_PATH})
      set(GRIDPACK_TEST_LDPATH "${GRIDPACK_TEST_LDPATH}:$ENV{LD_LIBRARY_PATH}")
    endif()
  endif()
  message(STATUS "GRIDPACK_TEST_LDPATH=${GRIDPACK_TEST_LDPATH}")
endif()

# -------------------------------------------------------------
# Installation
# -------------------------------------------------------------

install(FILES
  ${GRIDPACK_SRC_DIR}/../README.md
  LICENSE
  DESTINATION share/gridpack
)

# -------------------------------------------------------------
# Subdirectories
#
#  ****NOTE****
#  For any subdirectory that contributes a LIBRARY to the
#  GridPACK framework, modifications also need to be made to
#  the src/lib/GridPACK.cmake.in file so that installs will
#  work correctly
#
# -------------------------------------------------------------
add_subdirectory(utilities)
add_subdirectory(parallel)
add_subdirectory(environment)
add_subdirectory(configuration)
add_subdirectory(math)
add_subdirectory(component)
add_subdirectory(analysis)
add_subdirectory(partition)
add_subdirectory(factory)
add_subdirectory(mapper)
add_subdirectory(expression)
add_subdirectory(optimization)
add_subdirectory(applications/powerflow)
add_subdirectory(applications/dynamic_simulation_full_y)
add_subdirectory(applications/contingency_analysis)
add_subdirectory(applications/rtpr)
add_subdirectory(applications/state_estimation)
add_subdirectory(applications/kalman_ds)
add_subdirectory(applications/development/powerflow2)
add_subdirectory(applications/development/dynamic_simulation_reduced_y)
add_subdirectory(applications/development/dynamic_simulation_implicit)
add_subdirectory(applications/examples/hello_world)
add_subdirectory(applications/examples/resistor_grid)
add_subdirectory(applications/examples/powerflow)
add_subdirectory(applications/examples/contingency_analysis)
add_subdirectory(applications/modules/powerflow)
add_subdirectory(applications/modules/dynamic_simulation)
add_subdirectory(applications/modules/dynamic_simulation_full_y)
add_subdirectory(applications/modules/state_estimation)
add_subdirectory(applications/modules/kalman_ds)
add_subdirectory(applications/components/kds_matrix)
add_subdirectory(applications/components/ds_matrix)
add_subdirectory(applications/components/se_matrix)
add_subdirectory(applications/components/pf_matrix)
add_subdirectory(applications/components/y_matrix)
add_subdirectory(network)
add_subdirectory(stream)
add_subdirectory(parser)
add_subdirectory(export)
add_subdirectory(serial_io)
add_subdirectory(timer)
add_subdirectory(include)

add_subdirectory(lib)


# -------------------------------------------------------------
# CPack configuration
# -------------------------------------------------------------


INCLUDE(InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY 
  "An open source toolkit for developing power grid simulation applications for high performance computing architectures")
set(CPACK_PACKAGE_DESCRIPTION "GridPACK(tm) is a software framework consisting of a set of modules designed to simplify the development of programs that model the power grid and run on parallel, high performance computing platforms. The modules are available as a library and consist of components for setting up and distributing power grid networks, support for modeling the behavior of individual buses and branches in the network, converting the network models to the corresponding algebraic equations, and parallel routines for manipulating and solving large algebraic systems. Additional modules support input and output as well as basic profiling and error management. ")
set(CPACK_PACKAGE_VENDOR "")
set(CPACK_PACKAGE_CONTACT "william.perkins@pnnl.gov")
set(CPACK_PACKAGE_VERSION_MAJOR "${GridPACK_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${GridPACK_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${GridPACK_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION 
  "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_PACKAGE_FILE_NAME 
  "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME 
  "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION}")

# set(CPACK_PACKAGING_INSTALL_PREFIX "")
# set(CMAKE_INSTALL_PREFIX "/")
#set(CPACK_PACKAGE_INSTALL_DIRECTORY "/")
#set(CPACK_SET_DESTDIR FALSE)
set(CPACK_SET_DESTDIR true)

# Debian specific settings

set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})
set(CPACK_DEBIAN_PACKAGE_DEPENDS 
  "libga-dev (>=5.3), libpetsc3.6.2-dev (>= 3.6.2), libparmetis-dev (>=4.0.3), libboost-mpi-dev (>=1.58.0), libboost-all-dev (>=1.58.0)"
  )

# RPM specific settings

set(CPACK_RPM_PACKAGE_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})

INCLUDE(CPack)
