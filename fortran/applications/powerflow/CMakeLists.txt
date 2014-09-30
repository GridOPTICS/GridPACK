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
# Last Change: 2014-09-29 12:29:27 d3g096
# -------------------------------------------------------------

add_executable(powerflow.x
  pf_components.F90
  pf_factory.F90
  pf_main.F90
  $<TARGET_OBJECTS:gridpackf_network>
  $<TARGET_OBJECTS:gridpackf_factory>
  $<TARGET_OBJECTS:gridpackf_component>
  $<TARGET_OBJECTS:gridpackf_mapper>
  $<TARGET_OBJECTS:gridpackf_parallel>
  $<TARGET_OBJECTS:gridpackf_math>
  $<TARGET_OBJECTS:gridpackf_parser>
  $<TARGET_OBJECTS:gridpackf_serial_io>
  $<TARGET_OBJECTS:gridpackf_configuration>
)

target_link_libraries(powerflow.x
  ${MPI_Fortran_LIBRARIES}
  ${GRIDPACK_LIBS}
)

add_custom_target(powerflow_input

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/IEEE14.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS

  ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/IEEE14.raw
)

add_dependencies(powerflow.x
  powerflow_input
  gridpackf_network
  gridpackf_factory
  gridpackf_mapper
  gridpackf_component
  gridpackf_math
  gridpackf_parser
  gridpackf_serial_io
  gridpackf_configuration
  gridpackf_parallel
)
