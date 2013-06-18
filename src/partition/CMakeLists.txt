# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created June 17, 2013 by William A. Perkins
# Last Change: 2013-06-18 09:41:16 d3g096
# -------------------------------------------------------------


set(gridpack_partition_sources)

set(target_libraries gridpack_parallel ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})


# implementation independent partition sources
list(APPEND gridpack_partition_sources
  adjacency_list.cpp
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

# -------------------------------------------------------------
# implementation specific sources and flags
# -------------------------------------------------------------


# -------------------------------------------------------------
# gridpack partition library
# -------------------------------------------------------------
add_library(gridpack_partition
  ${gridpack_partition_sources}
)

# -------------------------------------------------------------
# gridpack partition header installation
# -------------------------------------------------------------

# -------------------------------------------------------------
# gridpack partition library installation
# -------------------------------------------------------------

# -------------------------------------------------------------
# partition test suite
# -------------------------------------------------------------
add_executable(partition_test 
  test/partition_test.cpp
)
target_link_libraries(partition_test gridpack_partition ${target_libraries})
gridpack_add_unit_test(partition partition_test)
