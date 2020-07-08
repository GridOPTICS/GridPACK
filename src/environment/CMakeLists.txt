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
# Created May  6, 2013 by William A. Perkins
# Last Change: 2020-05-22 12:12:23 d3g096
# -------------------------------------------------------------

add_library(gridpack_environment 
  environment.cpp
  no_print.cpp
  )

target_link_libraries(gridpack_environment
  gridpack_math
  )
gridpack_set_library_version(gridpack_environment)

set(target_libraries
   gridpack_environment
   gridpack_math
   ${GA_LIBRARIES}
   ${MPI_CXX_LIBRARIES}
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# -------------------------------------------------------------
# gridpack header installation
# -------------------------------------------------------------
install(FILES
  environment.hpp
  no_print.hpp
  DESTINATION include/gridpack/environment
)

# -------------------------------------------------------------
# gridpack environment library installation
# -------------------------------------------------------------
install(TARGETS 
  gridpack_environment
  DESTINATION lib
)

# -------------------------------------------------------------
# TEST: env_test
# A simple program to test the run time options capabilities
# of the environment module
# -------------------------------------------------------------
add_executable(env_test test/env_test.cpp)
target_link_libraries(env_test ${target_libraries})

gridpack_add_run_test(env_test env_test "-i input -o output")
