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
# Last Change: 2019-01-04 08:39:22 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_dynamic_simulation_module
    gridpack_powerflow_module
    gridpack_dsmatrix_components
    gridpack_pfmatrix_components
    gridpack_ymatrix_components
    gridpack_components
    gridpack_stream
    gridpack_partition
    gridpack_environment
    gridpack_math
    gridpack_configuration
    gridpack_timer
    gridpack_parallel
)

if (GOSS_DIR)
  set(target_libraries
      ${target_libraries}
      gridpack_goss
      ${GOSS_LIBRARY}
      ${APR_LIBRARY})
endif()

set (target_libraries
    ${target_libraries}
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_executable(ds.x
   ds_main.cpp
)

target_link_libraries(ds.x ${target_libraries})

add_custom_target(ds.x.input
 
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

add_dependencies(ds.x ds.x.input)

# -------------------------------------------------------------
# install as an example
# -------------------------------------------------------------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.install.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)

# install(FILES 
#   ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
#   ${CMAKE_CURRENT_BINARY_DIR}/input.xml
#   ${CMAKE_CURRENT_SOURCE_DIR}/IEEE_145bus_v23_PSLF.raw
#   ${CMAKE_CURRENT_SOURCE_DIR}/IEEE_145b_classical_model.dyr
#   ${CMAKE_CURRENT_SOURCE_DIR}/ds_main.cpp
#   DESTINATION share/gridpack/example/dynamic_simulation_reduced_y
# )

# install(TARGETS ds.x DESTINATION bin)
  

# -------------------------------------------------------------
# run application as test
# -------------------------------------------------------------
gridpack_add_run_test("dynamic_simulation" ds.x input.xml)

