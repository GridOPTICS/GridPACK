#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be
#     found
#     in the LICENSE file in the top level directory of this distribution.
#
# -*- mode: cmake -*-

# -------------------------------------------------------------
# file: CMakeLists.install.in
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created January 28, 2014 by William A. Perkins
# Last Change: 2014-09-30 13:52:17 d3g096
# -------------------------------------------------------------

cmake_minimum_required(VERSION 2.8.10)

if (NOT GRIDPACK_DIR)
  message(FATAL_ERROR
    "GRIDPACK_DIR variable must be set to the GridPACK installation directory")
endif()

include("${GRIDPACK_DIR}/lib/GridPACK.cmake")

project(GridPACKFortran)

enable_language(Fortran CXX)

gridpack_setup()

# -------------------------------------------------------------
# MPI compiler must be rechecked to get Fortran libraries
# -------------------------------------------------------------
message(STATUS "Checking MPI ...")
find_package (MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(AFTER ${MPI_Fortran_INCLUDE_PATH})

# -------------------------------------------------------------
# GridPACK definitions, include directories, and libraries
# -------------------------------------------------------------
add_definitions(${GRIDPACK_DEFINITIONS})
include_directories(BEFORE ${GRIDPACK_INCLUDE_DIRS})

# -------------------------------------------------------------
# Fortran compiler specific stuff 
# -------------------------------------------------------------

# GNU Fortran 4.4 (and perhaps other versions) does not understand the
# CLASS keyword as used in this project.  To get around this, the
# "CLASS" keyword is replaced with "TYPE"; on the CMake command line
# this option should appear:
#
#   -D CLASS_WORD:STRING="TYPE"
#
# It would be nice to just set this based on compiler version, but
# CMake has no way to find the Fortran compiler version.  They
# probably think Fortran is stupid ;->.
if (DEFINED CLASS_WORD)
  add_definitions(-DCLASS=${CLASS_WORD})
endif()

# -------------------------------------------------------------
# This project's settings
# -------------------------------------------------------------

# All Fortran modules are put in the same directory.
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Put the top source directory in the include directories
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

# Put the (Fortran) include directory in the include directories
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(parallel)
add_subdirectory(configuration)
add_subdirectory(math)
add_subdirectory(component)
add_subdirectory(network)
add_subdirectory(factory)
add_subdirectory(mapper)
add_subdirectory(parser)
add_subdirectory(serial_io)
add_subdirectory(applications/powerflow)
add_subdirectory(lib)
