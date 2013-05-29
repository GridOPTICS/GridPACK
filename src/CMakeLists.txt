# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May  3, 2013 by William A. Perkins
# Last Change: 2013-05-17 14:50:53 d3g096
# -------------------------------------------------------------

cmake_minimum_required (VERSION 2.6.4)
project (GridPACK)

# This is a C++ project exclusively
enable_language(CXX)

# where to look for special .cmake files
list (APPEND CMAKE_MODULE_PATH "${GridPACK_SOURCE_DIR}/cmake-modules")

# turn testing on
enable_testing()

# -------------------------------------------------------------
# MPI compiler
# -------------------------------------------------------------
message(STATUS "Checking MPI ...")
find_package (MPI REQUIRED)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

# -------------------------------------------------------------
# Boost
#
# Find the Boost library. It needs to use the same compiler as
# MPI. How can I check that?
# -------------------------------------------------------------
message(STATUS "Checking Boost ...")
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.52 COMPONENTS mpi REQUIRED)

# -------------------------------------------------------------
# PETSc
#
# We need a very recent version of PETSc. It needs to parallel and
# compiled with C++ complex support.
# -------------------------------------------------------------
message(STATUS "Checking PETSc ...")
find_package (PETSc COMPONENTS CXX REQUIRED)

if (NOT EXISTS "${PETSC_DIR}/${PETSC_ARCH}/conf/PETScConfig.cmake")
  message(FATAL_ERROR "PETSc found, but CMake configuration for PETSc installation not found?")
endif()

include("${PETSC_DIR}/${PETSC_ARCH}/conf/PETScConfig.cmake")

# checks 

if (NOT PETSC_HAVE_MPI)
  message(FATAL_ERROR "PETSc installation is not parallel (--with-mpi=1)")
endif()

if (NOT PETSC_USE_COMPLEX) 
  message (FATAL_ERROR "PETSc installation does not use complex (--with-scalar-type=complex)")
endif()

# This is probably redundant
if (NOT PETSC_CLANGUAGE_Cxx)
  message (FATAL_ERROR "PETSc installation does not use C++ (--with-clanguage=c++)")
endif()

# -------------------------------------------------------------
# Doxygen
# -------------------------------------------------------------
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
      )
endif(DOXYGEN_FOUND)

# -------------------------------------------------------------
# Options for all Programs
# -------------------------------------------------------------
include_directories(AFTER ${Boost_INCLUDE_DIRS})
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

# -------------------------------------------------------------
# Subdirectories
# -------------------------------------------------------------
add_subdirectory(parallel)
add_subdirectory(math)
add_subdirectory(network)
