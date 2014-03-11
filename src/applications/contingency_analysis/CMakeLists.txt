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
# Last Change: 2014-02-12 08:34:18 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_components
    gridpack_factory
    gridpack_partition
    gridpack_parallel
    gridpack_math
    gridpack_configuration
    gridpack_timer
    gridpack_pfmatrix_components
    gridpack_ymatrix_components
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${PETSC_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_executable(ca.x
   ca_driver.cpp
   ca_app.cpp
   ca_factory.cpp
   ca_components.cpp
   ca_main.cpp
)

target_link_libraries(ca.x ${target_libraries})

# Put some sample input in the binary directory so ca.x can run

add_custom_target(ca.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/SJin_3g9b.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/SJin_3g9b.raw
)
add_dependencies(ca.x ca.x.input)

# -------------------------------------------------------------
# TEST: pf_network_test
# -------------------------------------------------------------
#add_executable(pf_network_test test/pf_network_test.cpp)
#target_link_libraries(pf_network_test powerflow_app ${target_libraries})
#gridpack_add_unit_test(pf_network_test pf_network_test)

