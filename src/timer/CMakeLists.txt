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
# Last Change: 2016-10-03 13:44:18 d3g096
# -------------------------------------------------------------

add_library(gridpack_timer
  coarse_timer.cpp
  local_timer.cpp
)

# target_link_libraries(gridpack_timer
#     gridpack_parallel
#     ${GA_LIBRARIES}
#     ${Boost_LIBRARIES}
#     ${MPI_CXX_LIBRARIES}
# )

set(target_libraries
    gridpack_timer
    gridpack_parallel
    ${GA_LIBRARIES}
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(AFTER ${GA_INCLUDE_DIRS})

# -------------------------------------------------------------
# TEST: test_timer
# -------------------------------------------------------------

add_executable(test_timer test/test_timer.cpp)
target_link_libraries(test_timer ${target_libraries})

gridpack_add_unit_test(timer_test test_timer)

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  coarse_timer.hpp
  local_timer.hpp
  DESTINATION include/gridpack/timer
)

install(TARGETS 
  gridpack_timer
  DESTINATION lib
)
