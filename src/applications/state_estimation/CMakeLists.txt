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
# Last Change: 2017-12-08 08:40:35 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_state_estimation_module
    gridpack_sematrix_components
    gridpack_ymatrix_components
    gridpack_components
    gridpack_stream
    gridpack_partition
    gridpack_environment
    gridpack_math
    gridpack_configuration
    gridpack_timer
    gridpack_parallel
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    )

if (GOSS_DIR)
  set(target_libraries
      ${target_libraries}
      gridpack_goss
      ${GOSS_LIBRARY}
      ${APR_LIBRARY})
endif()

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_executable(stes.x
   se_main.cpp
)

target_link_libraries(stes.x ${target_libraries})

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_14.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/se/input_14.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_14.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/se/input_14.xml"
  )

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_118.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/se/input_118.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_118.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/se/input_118.xml"
  )

add_custom_target(stes.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/IEEE14.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/measurements/IEEE14_meas.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/IEEE118.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/measurements/IEEE118_meas.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_BINARY_DIR}/input_14.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE14.raw
  ${GRIDPACK_DATA_DIR}/measurements/IEEE14_meas.xml
  ${CMAKE_CURRENT_BINARY_DIR}/input_118.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE118.raw
  ${GRIDPACK_DATA_DIR}/measurements/IEEE118_meas.xml
)

add_dependencies(stes.x stes.x.input)

# -------------------------------------------------------------
# install as an example
# -------------------------------------------------------------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.install.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)

install(FILES 
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
  ${CMAKE_CURRENT_BINARY_DIR}/input_14.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE14.raw
  ${GRIDPACK_DATA_DIR}/measurements/IEEE14_meas.xml
  ${CMAKE_CURRENT_BINARY_DIR}/input_118.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE118.raw
  ${GRIDPACK_DATA_DIR}/measurements/IEEE118_meas.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/se_main.cpp
  DESTINATION share/gridpack/example/state_estimation
)

install(TARGETS stes.x DESTINATION bin)
  
# -------------------------------------------------------------
# run application as test
# -------------------------------------------------------------
gridpack_add_run_test("state_estimation" stes.x input_14.xml)
