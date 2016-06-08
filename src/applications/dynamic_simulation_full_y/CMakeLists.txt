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
# Last Change: 2016-06-08 12:20:53 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_dynamic_simulation_full_y_module
    gridpack_powerflow_module
    gridpack_dsmatrix_components
    gridpack_ymatrix_components
    gridpack_components
    gridpack_factory
    gridpack_partition
    gridpack_math
    gridpack_configuration
    gridpack_timer
    gridpack_parallel
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${PETSC_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_executable(dsf.x
   dsf_main.cpp
)

target_link_libraries(dsf.x ${target_libraries})

add_custom_target(dsf.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/IEEE_145bus_v23_PSLF.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/IEEE_145b_classical_model.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/IEEE_145bus_v23_PSLF.raw
  ${CMAKE_CURRENT_SOURCE_DIR}/IEEE_145b_classical_model.dyr
)

add_dependencies(dsf.x dsf.x.input)

# -------------------------------------------------------------
# run application as test
# -------------------------------------------------------------
gridpack_add_run_test("dynamic_simulation_full_y" dsf.x input.xml)

