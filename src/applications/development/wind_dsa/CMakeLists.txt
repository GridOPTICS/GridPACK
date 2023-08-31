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
# Last Change: 2023-08-31 09:05:17 d3g096
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
    gridpack_environment
    gridpack_math
    gridpack_configuration
    gridpack_timer
    gridpack_parallel
    gridpack_block_parsers
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

if (HELICS_DIR)
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
if (PARMETIS_FOUND)
  include_directories(BEFORE ${PARMETIS_INCLUDE_DIR})
endif()

add_executable(wind.x
  wind_main.cpp
  wind_driver.cpp
)

add_executable(wind2.x
  wind_main2.cpp
  wind_driver.cpp
  wind_driver2.cpp
)
target_link_libraries(wind.x ${target_libraries})
target_link_libraries(wind2.x ${target_libraries})


gridpack_set_lu_solver(
  "${CMAKE_CURRENT_SOURCE_DIR}/input.xml"
  "${CMAKE_CURRENT_BINARY_DIR}/input.xml"
)

gridpack_set_lu_solver(
  "${CMAKE_CURRENT_SOURCE_DIR}/input_tamu2000_dsf.xml"
  "${CMAKE_CURRENT_BINARY_DIR}/input_tamu2000_dsf.xml"
)

add_custom_target(wind.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/IEEE3G9B_V23.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/3g9b_faults.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/3g9b_classical.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/load.txt
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/wind.txt
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/texas2000_hr_4882_mod_33.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/texas2000_v23_fullsys_noac_wind_mod.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/tamu_2000_faults.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/texas_2k_scenfile_5_full.txt
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_BINARY_DIR}/input.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/IEEE3G9B_V23.raw
  ${CMAKE_CURRENT_SOURCE_DIR}/3g9b_faults.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/3g9b_classical.dyr
  ${CMAKE_CURRENT_SOURCE_DIR}/load.txt
  ${CMAKE_CURRENT_SOURCE_DIR}/wind.txt
  ${CMAKE_CURRENT_BINARY_DIR}/input_tamu2000_dsf.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/texas2000_hr_4882_mod_33.raw
  ${CMAKE_CURRENT_SOURCE_DIR}/texas2000_v23_fullsys_noac_wind_mod.dyr
  ${CMAKE_CURRENT_SOURCE_DIR}/tamu_2000_faults.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/texas_2k_scenfile_5_full.txt
)
add_dependencies(wind.x wind.x.input)
add_dependencies(wind2.x wind.x.input)

# -------------------------------------------------------------
# run applications as tests
# -------------------------------------------------------------
gridpack_add_run_test("wind_dsa" wind.x input.xml)
gridpack_add_run_test("wind_dsa_2" wind2.x input.xml)
gridpack_add_run_test("wind_dsa_tamu2000" wind.x input_tamu2000_dsf.xml)
gridpack_add_run_test("wind_dsa_2_tamu2000" wind2.x input_tamu2000_dsf.xml)

# -------------------------------------------------------------
# install as a sample application
# -------------------------------------------------------------

#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.install.in
#  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)

# install(FILES 
#   ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
#   ${CMAKE_CURRENT_SOURCE_DIR}/wind_main.cpp
#   ${CMAKE_CURRENT_SOURCE_DIR}/wind_driver.cpp
#   ${CMAKE_CURRENT_SOURCE_DIR}/wind_driver.hpp  
#   ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
#   ${CMAKE_CURRENT_SOURCE_DIR}/IEEE14.raw
#   ${CMAKE_CURRENT_SOURCE_DIR}/118.xml
#   ${CMAKE_CURRENT_SOURCE_DIR}/118.raw
#   DESTINATION share/gridpack/example/wind_dsa
# )
