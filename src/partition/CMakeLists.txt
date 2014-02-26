#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be
#     found
#     in the LICENSE file in the top level directory of this distribution.
#
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created June 17, 2013 by William A. Perkins
# Last Change: 2014-02-26 07:54:54 d3g096
# -------------------------------------------------------------


set(gridpack_partition_sources)

set(target_libraries 
  gridpack_parallel 
  gridpack_timer
  ${GA_LIBRARIES} 
  ${Boost_LIBRARIES} 
  ${MPI_CXX_LIBRARIES}
)


# implementation independent partition sources
list(APPEND gridpack_partition_sources
  adjacency_list.cpp
  graph_partitioner.cpp
  graph_partitioner_implementation.cpp
  simple_adjacency.cpp
)

# implementation independent headers
set(gridpack_partition_headers
  graph_partitioner_implementation.hpp
  adjacency_list.hpp
  simple_adjacency.hpp
  graph_partitioner.hpp
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(AFTER ${GA_INCLUDE_DIRS})

# -------------------------------------------------------------
# implementation specific sources and flags
# -------------------------------------------------------------
if (PARMETIS_FOUND)

   include_directories(BEFORE ${PARMETIS_INCLUDE_DIR})
   list(APPEND gridpack_partition_sources
     parmetis/graph_partitioner.cpp
     parmetis/parmetis_graph_wrapper.cpp
     parmetis/parmetis_graph_partitioner_impl.cpp
   )
   set(target_libraries 
     ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
     ${GA_LIBRARIES} 
     ${target_libraries}
   )

else()

    message(FATAL_ERROR "A partition implementation (e.g. ParMETIS) must be specified")

endif()

# -------------------------------------------------------------
# gridpack partition library
# -------------------------------------------------------------
add_library(gridpack_partition
  ${gridpack_partition_sources}
)

# -------------------------------------------------------------
# partition test suite
# -------------------------------------------------------------
add_executable(partition_test 
  test/partition_test.cpp
)
target_link_libraries(partition_test gridpack_partition ${target_libraries})
gridpack_add_unit_test(partition partition_test)



# -------------------------------------------------------------
# ParMETIS specific unit tests
# -------------------------------------------------------------
if (PARMETIS_FOUND) 
    add_executable(parmetis_test 
      test/parmetis_test.cpp
    )
    target_link_libraries(parmetis_test gridpack_partition ${target_libraries})
    gridpack_add_unit_test(parmetis parmetis_test)
endif (PARMETIS_FOUND)

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES
  ${gridpack_partition_headers}
  DESTINATION include/gridpack/partition
)

# -------------------------------------------------------------
# gridpack partition library installation
# -------------------------------------------------------------
install(TARGETS
  gridpack_partition
  DESTINATION lib
)


