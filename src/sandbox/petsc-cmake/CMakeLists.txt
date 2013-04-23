# -*- mode: cmake -*-

cmake_minimum_required (VERSION 2.6.3)
project (PETScTest)

list (APPEND CMAKE_MODULE_PATH "${PETScTest_SOURCE_DIR}/cmake-modules")

# This is a C++ project exclusively
enable_language(CXX)


# -------------------------------------------------------------
# MPI
# -------------------------------------------------------------
message(STATUS "Checking MPI ...")
find_package (MPI REQUIRED)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

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
# Boost
#
# Find the Boost library. It needs to use the same compiler as
# PETSc. How can I check that?
# -------------------------------------------------------------
message(STATUS "Checking Boost ...")
find_package(Boost 1.52 COMPONENTS mpi REQUIRED)

# -------------------------------------------------------------
# Options for all Programs
# -------------------------------------------------------------
include_directories (${PETSC_INCLUDES} ${Boost_INCLUDE_DIR})
add_definitions (${PETSC_DEFINITIONS})

# -------------------------------------------------------------
# Programs
# -------------------------------------------------------------
add_executable(vwrap vwrap.cpp)
target_link_libraries(vwrap ${PETSC_LIBRARIES} ${Boost_LIBRARIES})


