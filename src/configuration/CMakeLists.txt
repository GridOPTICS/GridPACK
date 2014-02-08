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
# Last Change: 2014-01-09 09:16:03 d3g096
# -------------------------------------------------------------

add_library(gridpack_configuration
  configuration.cpp
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(AFTER ${GA_INCLUDE_DIRS})

# -------------------------------------------------------------
# TEST: configuration_test
# -------------------------------------------------------------
# put unit test configuration file in build directory so it can be
# read there
add_custom_target(configuration_test_input
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/test/configuration_test.xml
  ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/test/configuration_test.xml
)

add_executable(configuration_test test/configuration_test.cpp)
target_link_libraries(configuration_test 
  gridpack_configuration
  ${GA_LIBRARIES}
  ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})
gridpack_add_unit_test(configuration_test configuration_test)

add_dependencies(configuration_test configuration_test_input)

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  configuration.hpp
  configurable.hpp
  DESTINATION include/gridpack/configuration
)

install(TARGETS 
  gridpack_configuration
  DESTINATION lib
)
