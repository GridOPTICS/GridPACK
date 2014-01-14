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
# Last Change: 2014-01-08 09:06:28 d3g096
# -------------------------------------------------------------

set(target_libraries 
  ${GA_LIBRARIES} 
  ${Boost_LIBRARIES} 
  ${MPI_CXX_LIBRARIES}
)

if (PARMETIS_FOUND)
  include_directories(BEFORE ${PARMETIS_INCLUDE_DIR})
  set(target_libraries 
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${target_libraries}
    )
endif()

set(target_libraries 
  gridpack_partition gridpack_components gridpack_parallel
  ${target_libraries} 
)

include_directories(AFTER ${GA_INCLUDE_DIRS})

#add_library(gridpack_network
#  test_network.cpp
#)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# -------------------------------------------------------------
# TEST: test_network
# -------------------------------------------------------------
add_executable(test_network test/test_network.cpp)
target_link_libraries(test_network ${target_libraries})
gridpack_add_unit_test(test_network test_network)

# -------------------------------------------------------------
# TEST: network_partition
# -------------------------------------------------------------
add_executable(network_partition test/network_partition.cpp)
target_link_libraries(network_partition ${target_libraries})
gridpack_add_unit_test(network_partition network_partition)

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  base_network.hpp
  DESTINATION include/gridpack/network
)

