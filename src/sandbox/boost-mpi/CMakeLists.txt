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
# Last Change: 2018-09-20 14:27:44 d3g096
# -------------------------------------------------------------

cmake_minimum_required (VERSION 2.6.4)
project (BoostMPI)

# This is a C++ project exclusively
enable_language(CXX)

# turn testing on
enable_testing()

# -------------------------------------------------------------
# MPI compiler
# -------------------------------------------------------------
message(STATUS "Checking MPI ...")
find_package (MPI REQUIRED)
include_directories(AFTER ${MPI_INCLUDE_PATH})

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
# MPI. How can I check that?
# -------------------------------------------------------------
message(STATUS "Checking Boost ...")
#set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.49 COMPONENTS mpi serialization random system REQUIRED)
include_directories(AFTER ${Boost_INCLUDE_DIRS})

# -------------------------------------------------------------
# Options for all Programs/Libraries
# -------------------------------------------------------------

set(BOOSTMPI_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

# -------------------------------------------------------------
# Subdirectories
# -------------------------------------------------------------
add_subdirectory(boost_test)
