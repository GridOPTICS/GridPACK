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
# Last Change: 2020-07-08 09:39:26 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_dynamic_simulation_full_y_module
    gridpack_powerflow_module
    gridpack_pfmatrix_components
    gridpack_ymatrix_components
    gridpack_components
    gridpack_partition
    gridpack_environment
    gridpack_math
    gridpack_stream
    gridpack_configuration
    gridpack_timer
    gridpack_analysis
    gridpack_parallel
    gridpack_environment
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

add_executable(rtpr.x
   rtpr_driver.cpp
   rtpr_main.cpp
)

target_link_libraries(rtpr.x ${target_libraries})

# Put some sample input in the binary directory so rtpr.x can run

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_14_classic.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/rtpr/input_14_classic.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_14_classic.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/rtpr/input_14_classic.xml"
  )

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_tamu500.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/rtpr/input_tamu500.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_tamu500.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/rtpr/input_tamu500.xml"
  )

add_custom_target(rtpr.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/IEEE14_ca_mod_rate6.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/IEEE14.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/contingencies/contingencies_14.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/TAMU_500_rmsmallgen_v23.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/tamu_500bus_detail.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS
  ${CMAKE_CURRENT_BINARY_DIR}/input_14_classic.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE14_ca_mod_rate6.raw
  ${GRIDPACK_DATA_DIR}/dyr/IEEE14.dyr
  ${CMAKE_CURRENT_BINARY_DIR}/input_tamu500.xml
  ${GRIDPACK_DATA_DIR}/raw/TAMU_500_rmsmallgen_v23.raw
)
add_dependencies(rtpr.x rtpr.x.input)

# -------------------------------------------------------------
# install as an example
# -------------------------------------------------------------
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.install.in
#  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)

#install(FILES 
#  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
#  ${CMAKE_CURRENT_BINARY_DIR}/input_14_classic.xml
#  ${GRIDPACK_DATA_DIR}/raw/IEEE14_ca_mod_rate6.raw
#  ${GRIDPACK_DATA_DIR}/dyr/IEEE14.dyr
#  ${GRIDPACK_DATA_DIR}/contingencies/contingencies_14.xml
#  ${CMAKE_CURRENT_BINARY_DIR}/input_118.xml
#  ${GRIDPACK_DATA_DIR}/raw/IEEE118.raw
#  ${GRIDPACK_DATA_DIR}/contingencies/contingencies_118.xml
#  ${CMAKE_CURRENT_SOURCE_DIR}/rtpr_driver.cpp
#  ${CMAKE_CURRENT_SOURCE_DIR}/rtpr_driver.hpp
#  ${CMAKE_CURRENT_SOURCE_DIR}/rtpr_main.cpp
#  DESTINATION share/gridpack/example/rtpr
#)

install(TARGETS rtpr.x DESTINATION bin)

# -------------------------------------------------------------
# run application as test
# -------------------------------------------------------------
set(TIMEOUT 120.0)
gridpack_add_run_test("real_time_path_rating" rtpr.x input_14_classic.xml)

