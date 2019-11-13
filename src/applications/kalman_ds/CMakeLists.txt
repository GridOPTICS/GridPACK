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
# Last Change: 2018-07-17 13:41:33 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_kalmands_module
    gridpack_powerflow_module
    gridpack_kdsmatrix_components
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

include_directories(BEFORE
 ${CMAKE_CURRENT_SOURCE_DIR}/../modules/kalman_ds)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_executable(kds.x
   kds_main.cpp
)

target_link_libraries(kds.x ${target_libraries})


add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_145.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/kalman/input_145.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_145.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/kalman/input_145.xml"
  )

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_14.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/kalman/input_14.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_14.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/kalman/input_14.xml"
  )

add_custom_target(kds.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/IEEE145.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/IEEE145_classicGen.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/kalman/IEEE145_Kalman_input_ang.csv
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/kalman/IEEE145_Kalman_input_mag.csv
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/IEEE14_kds.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/IEEE14_classicGen.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/kalman/IEEE14_Kalman_input_ang.csv
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/kalman/IEEE14_Kalman_input_mag.csv
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_BINARY_DIR}/input_145.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE145.raw
  ${GRIDPACK_DATA_DIR}/dyr/IEEE145_classicGen.dyr
  ${GRIDPACK_DATA_DIR}/kalman/IEEE145_Kalman_input_ang.csv
  ${GRIDPACK_DATA_DIR}/kalman/IEEE145_Kalman_input_mag.csv
  ${CMAKE_CURRENT_BINARY_DIR}/input_14.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE14_kds.raw
  ${GRIDPACK_DATA_DIR}/dyr/IEEE14_classicGen.dyr
  ${GRIDPACK_DATA_DIR}/kalman/IEEE14_Kalman_input_ang.csv
  ${GRIDPACK_DATA_DIR}/kalman/IEEE14_Kalman_input_mag.csv
)

add_dependencies(kds.x kds.x.input)

# -------------------------------------------------------------
# install as an example
# -------------------------------------------------------------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.install.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)

install(FILES 
  ${CMAKE_CURRENT_SOURCE_DIR}/kds_main.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
  ${CMAKE_CURRENT_BINARY_DIR}/input_145.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE145.raw
  ${GRIDPACK_DATA_DIR}/dyr/IEEE145_classicGen.dyr
  ${GRIDPACK_DATA_DIR}/kalman/IEEE145_Kalman_input_ang.csv
  ${GRIDPACK_DATA_DIR}/kalman/IEEE145_Kalman_input_mag.csv
  ${CMAKE_CURRENT_BINARY_DIR}/input_14.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE14_kds.raw
  ${GRIDPACK_DATA_DIR}/dyr/IEEE14_classicGen.dyr
  ${GRIDPACK_DATA_DIR}/kalman/IEEE14_Kalman_input_ang.csv
  ${GRIDPACK_DATA_DIR}/kalman/IEEE14_Kalman_input_mag.csv
  DESTINATION share/gridpack/example/kalman_ds
)

install(TARGETS kds.x DESTINATION bin)
  


# -------------------------------------------------------------
# run application as test
# -------------------------------------------------------------
set(TIMEOUT 120.0)
gridpack_add_run_test("kalman_ds" kds.x input_14.xml)

