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
# Last Change: 2014-09-29 12:31:29 d3g096
# -------------------------------------------------------------

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(BEFORE ${GRIDPACK_INCLUDE_DIRS})

add_library(gridpackf_serial_io OBJECT
  serial_io_c.cpp
  serial_io_f.F90
)

add_dependencies(gridpackf_serial_io
  gridpackf_network
)

add_executable(fortran_serial_io_test
  serial_io_test.F90
  serial_io_component.F90
  serial_io_factory.F90
  $<TARGET_OBJECTS:gridpackf_serial_io>
  $<TARGET_OBJECTS:gridpackf_network>
  $<TARGET_OBJECTS:gridpackf_parallel>
  $<TARGET_OBJECTS:gridpackf_component>
  $<TARGET_OBJECTS:gridpackf_factory>
)

target_link_libraries(fortran_serial_io_test
  ${MPI_Fortran_LIBRARIES}
  ${GRIDPACK_LIBS}
)

