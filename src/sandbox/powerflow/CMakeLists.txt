# -*- mode: cmake -*-
#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be
#     found
#     in the LICENSE file in the top level directory of this distribution.
#
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May  6, 2013 by William A. Perkins
# Last Change: 2014-01-10 14:25:18 d3g096
# -------------------------------------------------------------
cmake_minimum_required(VERSION 2.6.4)

if (NOT GRIDPACK_DIR)
  message(FATAL_ERROR "You must specify a GridPACK installation with -D GRIDPACK_DIR=path")
endif()

set(GRIDPACK_CONFIG "${GRIDPACK_DIR}/lib/GridPACK.cmake")

if (NOT EXISTS ${GRIDPACK_CONFIG})
  message(FATAL_ERROR "GridPACK installation problem, cannot find ${GRIDPACK_CONFIG}")
endif()

include(${GRIDPACK_CONFIG})

project(PowerFlow)

enable_language(CXX)

gridpack_setup()

add_definitions(${GRIDPACK_DEFINITIONS})
include_directories(BEFORE ${GRIDPACK_INCLUDE_DIRS})

add_executable(powerflow.x
  pf_main.cpp
  pf_app.cpp
  pf_components.cpp
  pf_factory.cpp
)
target_link_libraries(powerflow.x ${GRIDPACK_LIBS})

add_executable(powerflow2
  pf_app2.cpp
  pf_components.cpp
  pf_factory.cpp
)
target_link_libraries(powerflow2 ${GRIDPACK_LIBS})

add_custom_target(powerflow.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/118.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/118.raw
)
add_dependencies(powerflow.x powerflow.input)
add_dependencies(powerflow2 powerflow.input)

