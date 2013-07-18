# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May  6, 2013 by William A. Perkins
# Last Change: 2013-07-18 09:00:33 d3g096
# -------------------------------------------------------------

set(target_libraries 
  gridpack_components gridpack_parallel
  ${GA_LIBRARIES} 
  ${Boost_LIBRARIES} 
  ${MPI_CXX_LIBRARIES}
)

include_directories(AFTER ${GA_INCLUDE_DIRS})

add_library(gridpack_network
  test_network.cpp
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# -------------------------------------------------------------
# TEST: greetings
# A simple MPI program
# -------------------------------------------------------------
add_executable(test_network test/test_network.cpp)
target_link_libraries(test_network gridpack_network ${target_libraries})

add_test(test_network_serial test_network)


add_executable(network_partition test/network_partition.cpp)
target_link_libraries(network_partition gridpack_network ${target_libraries})

gridpack_add_unit_test(network_partition network_partition)

