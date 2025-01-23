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
# Last Change: 2024-11-18 13:29:01 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_emt_module
    gridpack_powerflow_module
    gridpack_pfmatrix_components
    gridpack_ymatrix_components
    gridpack_components
    gridpack_partition
    gridpack_math
    gridpack_configuration
    gridpack_timer
    gridpack_parallel
    gridpack_block_parsers
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

include_directories(BEFORE
 ${CMAKE_CURRENT_SOURCE_DIR}/../modules/emt/model_classes)
include_directories(BEFORE
 ${CMAKE_CURRENT_SOURCE_DIR}/../modules/emt/base_classes)
include_directories(BEFORE
 ${CMAKE_CURRENT_SOURCE_DIR}/../modules/emt)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_executable(emt.x
   emt_app.cpp
)

target_link_libraries(emt.x ${target_libraries})

add_custom_target(emt.x.input

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/petscoptions/emt/.petscrc
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/emt/9-bus/input_9bus2.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/emt/9-bus/case9_PSCAD.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/emt/9-bus/case9_GAST.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/emt/39-bus/input_39bus.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/emt/39-bus/IEEE39bus.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/emt/39-bus/IEEE39bus.dyr
  ${CMAKE_CURRENT_BINARY_DIR}


  DEPENDS
  ${CMAKE_CURRENT_BINARY_DIR}/input_9bus2.xml
  ${CMAKE_CURRENT_BINARY_DIR}/input_39bus.xml
  ${GRIDPACK_DATA_DIR}/emt/9-bus/input_9bus2.xml
  ${GRIDPACK_DATA_DIR}/emt/39-bus/input_39bus.xml
  ${CMAKE_CURRENT_BINARY_DIR}/case9_PSCAD.raw
  ${CMAKE_CURRENT_BINARY_DIR}/IEEE39bus.raw
  ${GRIDPACK_DATA_DIR}/emt/9-bus/case9_PSCAD.raw
  ${GRIDPACK_DATA_DIR}/emt/39-bus/IEEE39bus.raw
  ${CMAKE_CURRENT_BINARY_DIR}/case9_GAST.dyr
  ${CMAKE_CURRENT_BINARY_DIR}/IEEE39bus.dyr
  ${GRIDPACK_DATA_DIR}/emt/9-bus/case9_GAST.dyr
  ${GRIDPACK_DATA_DIR}/emt/39-bus/IEEE39bus.dyr
  ${CMAKE_CURRENT_BINARY_DIR}/.petscrc
)

add_dependencies(emt.x emt.x.input)

install(FILES 
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
  ${CMAKE_CURRENT_BINARY_DIR}/input_9bus2.xml
  ${CMAKE_CURRENT_BINARY_DIR}/input_39bus.xml
  ${GRIDPACK_DATA_DIR}/emt/9-bus/input_9bus2.xml
  ${GRIDPACK_DATA_DIR}/emt/39-bus/input_39bus.xml
  ${CMAKE_CURRENT_BINARY_DIR}/case9_PSCAD.raw
  ${CMAKE_CURRENT_BINARY_DIR}/IEEE39bus.raw
  ${GRIDPACK_DATA_DIR}/emt/9-bus/case9_PSCAD.raw
  ${GRIDPACK_DATA_DIR}/emt/39-bus/IEEE39bus.raw
  ${CMAKE_CURRENT_BINARY_DIR}/case9_GAST.dyr
  ${CMAKE_CURRENT_BINARY_DIR}/IEEE39bus.dyr
  ${GRIDPACK_DATA_DIR}/emt/9-bus/case9_GAST.dyr
  ${GRIDPACK_DATA_DIR}/emt/39-bus/IEEE39bus.dyr
  ${GRIDPACK_DATA_DIR}/petscoptions/emt/.petscrc
  ${CMAKE_CURRENT_SOURCE_DIR}/emt_app.cpp
  DESTINATION share/gridpack/example/emt
)

install(TARGETS emt.x DESTINATION bin)

# -------------------------------------------------------------
# run application as test
# -------------------------------------------------------------
gridpack_add_run_test("EMT" emt.x input_9bus2.xml)

