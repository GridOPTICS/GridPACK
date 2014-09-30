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
# Created May 12, 2014 by William A. Perkins
# Last Change: 2014-08-22 08:08:26 d3g096
# -------------------------------------------------------------

add_library(gridpackf_configuration OBJECT
  configurationc.cpp
  configurationf.F90
)

add_dependencies(gridpackf_configuration
  gridpackf_parallel
)

add_executable(configuration_test
  configuration_test.F90
  $<TARGET_OBJECTS:gridpackf_configuration>
  $<TARGET_OBJECTS:gridpackf_parallel>
)
target_link_libraries(configuration_test
  ${GRIDPACK_LIBS}
)

add_custom_target(configuration_test_input

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/gridpack.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS

  ${CMAKE_CURRENT_SOURCE_DIR}/gridpack.xml
)
add_dependencies(configuration_test configuration_test_input)
