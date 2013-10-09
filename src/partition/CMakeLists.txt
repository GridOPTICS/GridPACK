# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created June 17, 2013 by William A. Perkins
# Last Change: 2013-10-09 14:59:20 d3g096
# -------------------------------------------------------------


set(gridpack_partition_sources)

set(target_libraries 
  gridpack_parallel 
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
# gridpack partition header installation
# -------------------------------------------------------------
install(FILES
  adjacency_list.hpp
  graph_partitioner.hpp
  graph_partitioner_implementation.hpp
  DESTINATION include/gridpack/partition
)

# -------------------------------------------------------------
# gridpack partition library installation
# -------------------------------------------------------------
install(TARGETS
  gridpack_partition
  DESTINATION lib
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