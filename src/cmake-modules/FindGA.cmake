# -*- mode: cmake -*-
# GridPACK Note:  This file is base on that found 18June2013 here:
# https://github.com/OP2/OP2-Common/blob/master/cmake/modules/FindParMETIS.cmake
#
# - Try to find GA libraries
# Once done this will define
#
#  GA_FOUND        - system has GA
#  GA_INCLUDE_DIRS - include directories for GA
#  GA_LIBRARIES    - libraries for GA
#
# Variables used by this module. They can change the default behaviour and
# need to be set before calling find_package:
#
#  GA_DIR          - Prefix directory of the GA installation
#  GA_INCLUDE_DIR  - Include directory of the GA installation
#                          (set only if different from ${GA_DIR}/include)
#  GA_LIB_DIR      - Library directory of the GA installation
#                          (set only if different from ${GA_DIR}/lib)
#  GA_TEST_RUNS    - Skip tests building and running a test
#                          executable linked against GA libraries
#  GA_LIB_SUFFIX   - Also search for non-standard library names with the
#                          given suffix appended

#=============================================================================
# Copyright (C) 2010-2012 Garth N. Wells, Anders Logg, Johannes Ring
# and Florian Rathgeber. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

# GA reqires MPI
find_package(MPI QUIET REQUIRED)

# Stupid RHEL puts GA headers in with the MPI headers

find_path(GA_INCLUDE_DIR ga++.h
  HINTS ${GA_INCLUDE_DIR} ${GA_DIR} ${MPI_CXX_INCLUDE_PATH}
  ENV GA_INCLUDE_DIR 
  ENV MPI_INCLUDE_DIR
  ENV GA_DIR
  PATH_SUFFIXES include ""
  DOC "Directory where the GA header files are located"
)
message(STATUS "GA_INCLUDE_DIR: ${GA_INCLUDE_DIR}")

# This is a hack for stupid RHEL rpms to look for the GA libraries in
# amongst the MPI libraries

if(MPI_CXX_LIBRARIES)
  list(GET MPI_CXX_LIBRARIES 0 lib1)
  get_filename_component(MPI_LIB_DIR_MAYBE ${lib1} DIRECTORY)
elseif(MPI_LIBRARY)
  get_filename_component(MPI_LIB_DIR_MAYBE ${MPI_LIBRARY} DIRECTORY)
endif()

# message(STATUS "MPI_LIB_DIR_MAYBE: ${MPI_LIB_DIR_MAYBE}")

find_library(GA_LIBRARY
  NAMES ga GA${GA_LIB_SUFFIX}
  HINTS ${GA_DIR} ${MPI_LIB_DIR_MAYBE}
  ENV GA_LIB_DIR ENV GA_DIR 
  ENV MPI_LIB
  PATH_SUFFIXES lib ""
  DOC "Directory where the GA library is located"
)
message(STATUS "GA_LIBRARY: ${GA_LIBRARY}")

find_library(GA_CXX_LIBRARY
  NAMES ga++ GA${GA_LIB_SUFFIX}
  HINTS ${GA_LIB_DIR} ${GA_DIR} ${MPI_LIB_DIR_MAYBE}
  ENV GA_LIB_DIR ENV GA_DIR
  PATH_SUFFIXES lib
  DOC "Directory where the GA library is located"
)
message(STATUS "GA_CXX_LIBRARY: ${GA_CXX_LIBRARY}")

find_library(ARMCI_LIBRARY
  NAMES armci GA${ARMCI_LIB_SUFFIX}
  HINTS ${GA_LIB_DIR} ${GA_DIR} ${MPI_LIB_DIR_MAYBE}
  ENV GA_LIB_DIR ENV GA_DIR
  PATH_SUFFIXES lib
  DOC "Directory where the GA library is located"
)
message(STATUS "ARMCI_LIBRARY: ${ARMCI_LIBRARY}")

# Get GA version
if(NOT GA_VERSION_STRING AND GA_INCLUDE_DIR AND EXISTS "${GA_INCLUDE_DIR}/ga.h")
  set(version_pattern "^#define[\t ]+GA_(MAJOR|MINOR)_VERSION[\t ]+([0-9\\.]+)$")
  file(STRINGS "${GA_INCLUDE_DIR}/gacommon.h" GA_version REGEX ${version_pattern})

  foreach(match ${GA_version})
    if(GA_VERSION_STRING)
      set(GA_VERSION_STRING "${GA_VERSION_STRING}.")
    endif()
    string(REGEX REPLACE ${version_pattern} "${GA_VERSION_STRING}\\2" GA_VERSION_STRING ${match})
    set(GA_VERSION_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
  endforeach()
  unset(GA_version)
  unset(version_pattern)
endif()

# Try compiling and running test program
if (GA_INCLUDE_DIR AND GA_LIBRARY AND ARMCI_LIBRARY)

  # Set flags for building test program
  set(CMAKE_REQUIRED_INCLUDES ${GA_INCLUDE_DIR} ${MPI_INCLUDE_PATH})
  if (NOT MPI_LIBRARY)
    set(CMAKE_REQUIRED_LIBRARIES 
      ${GA_CXX_LIBRARY} ${GA_LIBRARY} ${ARMCI_LIBRARY} ${GA_EXTRA_LIBS}
    )
  else()
    set(CMAKE_REQUIRED_LIBRARIES 
      ${GA_CXX_LIBRARY} ${GA_LIBRARY} ${ARMCI_LIBRARY} ${GA_EXTRA_LIBS} ${MPI_LIBRARIES}
    )
  endif()

# Build and run test program, maybe

set(ga_test_src "
#include <mpi.h>
#include <ga++.h>

int main(int argc, char **argv)
{
  // FIXME: Find a simple but sensible test for GA

  // Initialise MPI
  MPI_Init(&argc, &argv);

  // Initialize GA
  GA_Initialize();

  // Terminate GA
  GA_Terminate();

  // Finalize MPI
  MPI_Finalize();

  return 0;
}
")

include(CheckCXXSourceRuns)
include(CheckCXXSourceCompiles)
if (USE_PROGRESS_RANKS OR CHECK_COMPILATION_ONLY) 
  check_cxx_source_compiles("${ga_test_src}" GA_TEST_RUNS)
else()
  check_cxx_source_runs("${ga_test_src}" GA_TEST_RUNS)
endif()
endif()

unset(ga_test_src)

# Standard package handling
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.2)
  find_package_handle_standard_args(GA
    REQUIRED_VARS GA_LIBRARY GA_INCLUDE_DIR GA_TEST_RUNS
    VERSION_VAR GA_VERSION_STRING)
else()
  find_package_handle_standard_args(GA
    REQUIRED_VARS GA_LIBRARY GA_INCLUDE_DIR GA_TEST_RUNS)
endif()

if(GA_FOUND)
  set(GA_LIBRARIES ${GA_CXX_LIBRARY} ${GA_LIBRARY} ${ARMCI_LIBRARY} ${GA_EXTRA_LIBS})
  set(GA_INCLUDE_DIRS ${GA_INCLUDE_DIR})
endif()

mark_as_advanced(GA_INCLUDE_DIR GA_LIBRARY ARMCI_LIBRARY)
