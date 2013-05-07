# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May  3, 2013 by William A. Perkins
# Last Change: 2013-05-07 13:18:45 d3g096
# -------------------------------------------------------------

cmake_minimum_required (VERSION 2.6.3)
project (GridPACK)

# This is a C++ project exclusively
enable_language(CXX)

# where to look for special .cmake files
list (APPEND CMAKE_MODULE_PATH "${GridPACK_SOURCE_DIR}/cmake-modules")

# turn testingo on
enable_testing()

# -------------------------------------------------------------
# MPI compiler
# -------------------------------------------------------------
message(STATUS "Checking MPI ...")
find_package (MPI REQUIRED)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# -------------------------------------------------------------
# Boost
#
# Find the Boost library. It needs to use the same compiler as
# MPI. How can I check that?
# -------------------------------------------------------------
message(STATUS "Checking Boost ...")
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.52 COMPONENTS mpi REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# -------------------------------------------------------------
# Subdirectories
# -------------------------------------------------------------
add_subdirectory(parallel)
