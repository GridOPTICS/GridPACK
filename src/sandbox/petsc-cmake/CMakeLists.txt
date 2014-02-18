# -*- mode: cmake -*-

cmake_minimum_required (VERSION 2.6.4)
project (PETScTest)

list (APPEND CMAKE_MODULE_PATH "${PETScTest_SOURCE_DIR}/cmake-modules")

# This is a C++ project exclusively
enable_language(CXX)

# Do some tests
enable_testing()


# -------------------------------------------------------------
# MPI
# -------------------------------------------------------------
message(STATUS "Checking MPI ...")
find_package (MPI REQUIRED)

# This is here for cmake version 2.6 
if (NOT MPI_CXX_LIBRARIES)
 set(MPI_CXX_LIBRARIES 
   ${MPI_LIBRARY} 
   ${MPI_EXTRA_LIBRARY}
 )
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
# Options for all Programs
# -------------------------------------------------------------
include_directories (${PETSC_INCLUDES} ${Boost_INCLUDE_DIR})
add_definitions (${PETSC_DEFINITIONS})

# -------------------------------------------------------------
# Programs
# -------------------------------------------------------------
add_executable(vwrap vwrap.cpp)
target_link_libraries(vwrap ${PETSC_LIBRARIES} ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})

add_executable(vecdist vecdist.cpp)
target_link_libraries(vecdist ${PETSC_LIBRARIES} ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})

add_executable(exception exception.cpp)
target_link_libraries(exception ${PETSC_LIBRARIES} ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})


add_executable(greetings greetings.cpp)
target_link_libraries(greetings ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})
add_test(NAME greetings COMMAND $<TARGET_FILE:greetings>)
set_tests_properties(greetings
  PROPERTIES 
  PASS_REGULAR_EXPRESSION "I am process 0 of [1-9].*$"
)

