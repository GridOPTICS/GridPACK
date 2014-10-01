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
# Created September 9, 2014 by Bruce Palmer
# Last Change: 2014-09-29 15:03:31 d3g096
# -------------------------------------------------------------

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(BEFORE ${GRIDPACK_INCLUDE_DIRS})

add_library(gridpackf_parser OBJECT
  parser_c.cpp
  parser_f.F90
)

add_dependencies(gridpackf_parser
  gridpackf_network
)

add_executable(fortran_parser_test
  parser_test.F90
  parser_component.F90
  $<TARGET_OBJECTS:gridpackf_parser>
  $<TARGET_OBJECTS:gridpackf_network>
  $<TARGET_OBJECTS:gridpackf_parallel>
  $<TARGET_OBJECTS:gridpackf_component>
)

target_link_libraries(fortran_parser_test
  ${MPI_Fortran_LIBRARIES}
  ${GRIDPACK_LIBS}
)

add_custom_target(parser_input

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/parser_data.raw
  ${CMAKE_CURRENT_BINARY_DIR}

)

add_dependencies(fortran_parser_test
  parser_input
)
