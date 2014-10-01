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
# Created August 5, 2014 by Bruce Palmer
# Last Change: 2014-09-29 12:31:07 d3g096
# -------------------------------------------------------------

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(BEFORE ${GRIDPACK_INCLUDE_DIRS})

add_library(gridpackf_network OBJECT
  network_c.cpp
  network_f.F90
)

add_dependencies(gridpackf_network
  gridpackf_component
)

add_executable(fortran_network_test
  network_test.F90
  test_component.F90
  $<TARGET_OBJECTS:gridpackf_network>
  $<TARGET_OBJECTS:gridpackf_parallel>
  $<TARGET_OBJECTS:gridpackf_component>
)

target_link_libraries(fortran_network_test
  ${MPI_Fortran_LIBRARIES}
  ${GRIDPACK_LIBS}
)

