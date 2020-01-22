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
# Last Change: 2017-12-08 09:39:48 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_dynamic_simulation_full_y_module
    gridpack_powerflow_module
    gridpack_pfmatrix_components
    gridpack_dsmatrix_components
    gridpack_ymatrix_components
    gridpack_components
    gridpack_stream
    gridpack_partition
    gridpack_math
    gridpack_configuration
    gridpack_timer
    gridpack_parallel
    ${PETSC_LIBRARIES}
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

if (GOSS_DIR)
  set(target_libraries
      ${target_libraries}
      gridpack_goss
      ${GOSS_LIBRARY}
      ${APR_LIBRARY})
endif()

if (HELICS_INSTALL_DIR)
  set(target_libraries
      ${target_libraries}
      ${JSON_LIBRARY}
      ${ZEROMQ_LIBRARY}
      ${SODIUM_LIBRARY}
      ${HELICS_LIBRARY})
endif()

include_directories(BEFORE
 ${CMAKE_CURRENT_SOURCE_DIR}/../modules/dynamic_simulation_full_y/model_classes)
include_directories(BEFORE
 ${CMAKE_CURRENT_SOURCE_DIR}/../modules/dynamic_simulation_full_y/base_classes)
include_directories(BEFORE
 ${CMAKE_CURRENT_SOURCE_DIR}/../modules/dynamic_simulation_full_y)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_executable(dsf.x
   dsf_main.cpp
)

target_link_libraries(dsf.x ${target_libraries})


add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_145.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/ds/input_145.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_145.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/ds/input_145.xml"
  )

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_9b3g.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/ds/input_9b3g.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_9b3g.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/ds/input_9b3g.xml"
  )

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_300_cmpld.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/ds/input_300_cmpld.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_300_cmpld.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/ds/input_300_cmpld.xml"
  )

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_3000.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/ds/input_3000.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_3000.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/ds/input_3000.xml"
  )

add_custom_target(dsf.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/IEEE_145bus_v23_PSLF.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/IEEE_145b_classical_model.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/9b3g.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/9b3g.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/300bus_v23_no0imp_pslf.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/300bus_detail_model_cmpld_combine.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/bus3000_gen_no0imp_v23_pslf.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/classical_model_3000bus.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_BINARY_DIR}/input_145.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE_145bus_v23_PSLF.raw
  ${GRIDPACK_DATA_DIR}/dyr/IEEE_145b_classical_model.dyr
  ${CMAKE_CURRENT_BINARY_DIR}/input_9b3g.xml
  ${GRIDPACK_DATA_DIR}/raw/9b3g.raw
  ${GRIDPACK_DATA_DIR}/dyr/9b3g.dyr
  ${CMAKE_CURRENT_BINARY_DIR}/input_300_cmpld.xml
  ${GRIDPACK_DATA_DIR}/raw/300bus_v23_no0imp_pslf.raw
  ${GRIDPACK_DATA_DIR}/dyr/300bus_detail_model_cmpld_combine.dyr
  ${CMAKE_CURRENT_BINARY_DIR}/input_3000.xml
  ${GRIDPACK_DATA_DIR}/raw/bus3000_gen_no0imp_v23_pslf.raw
  ${GRIDPACK_DATA_DIR}/dyr/classical_model_3000bus.dyr
)

add_dependencies(dsf.x dsf.x.input)

# -------------------------------------------------------------
# install as an example
# -------------------------------------------------------------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.install.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)

install(FILES 
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
  ${CMAKE_CURRENT_BINARY_DIR}/input_145.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE_145bus_v23_PSLF.raw
  ${GRIDPACK_DATA_DIR}/dyr/IEEE_145b_classical_model.dyr
  ${CMAKE_CURRENT_BINARY_DIR}/input_9b3g.xml
  ${GRIDPACK_DATA_DIR}/raw/9b3g.raw
  ${GRIDPACK_DATA_DIR}/dyr/9b3g.dyr
  ${CMAKE_CURRENT_BINARY_DIR}/input_300_cmpld.xml
  ${GRIDPACK_DATA_DIR}/raw/300bus_v23_no0imp_pslf.raw
  ${GRIDPACK_DATA_DIR}/dyr/300bus_detail_model_cmpld_combine.dyr
  ${CMAKE_CURRENT_BINARY_DIR}/input_3000.xml
  ${GRIDPACK_DATA_DIR}/raw/bus3000_gen_no0imp_v23_pslf.raw
  ${GRIDPACK_DATA_DIR}/dyr/classical_model_3000bus.dyr
  ${CMAKE_CURRENT_SOURCE_DIR}/dsf_main.cpp
  DESTINATION share/gridpack/example/dynamic_simulation_full_y
)

install(TARGETS dsf.x DESTINATION bin)
  


# -------------------------------------------------------------
# run application as test
# -------------------------------------------------------------
gridpack_add_run_test("dynamic_simulation_full_y" dsf.x input_145.xml)

