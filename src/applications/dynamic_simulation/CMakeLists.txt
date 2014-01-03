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
# Last Change: 2013-08-09 11:41:42 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_components
    gridpack_factory
    gridpack_partition
    gridpack_parallel
    gridpack_math
    gridpack_configuration
    gridpack_timer
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${PETSC_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_executable(ds.x
   ds_app.cpp
   ds_components.cpp
   ds_factory.cpp
   ds_main.cpp
)

target_link_libraries(ds.x ${target_libraries})


# -------------------------------------------------------------
# TEST: pf_network_test
# -------------------------------------------------------------
#add_executable(pf_network_test test/pf_network_test.cpp)
#target_link_libraries(pf_network_test powerflow_app ${target_libraries})
#gridpack_add_unit_test(pf_network_test pf_network_test)

