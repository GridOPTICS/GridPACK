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
# Last Change: 2020-05-22 09:41:25 d3g096
# -------------------------------------------------------------

set(target_libraries
   gridpack_parallel
   gridpack_math
   gridpack_timer
   gridpack_environment
   gridpack_math
   ${Boost_LIBRARIES} 
   ${GA_LIBRARIES} 
   ${MPI_CXX_LIBRARIES}
)
	
add_library(gridpack_parallel 
  communicator.cpp
  distributed.cpp
  index_hash.cpp
  random.cpp
  )
add_dependencies(gridpack_parallel external_build)
gridpack_set_library_version(gridpack_parallel)

# For some reason this is necessary for Mac OS X, but messes things up on Linux
target_link_libraries(gridpack_parallel
  ${Boost_LIBRARIES} 
  ${GA_LIBRARIES} 
  ${MPI_CXX_LIBRARIES}
  )

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# -------------------------------------------------------------
# gridpack header installation
# -------------------------------------------------------------
install(FILES
  communicator.hpp
  parallel.hpp
  distributed.hpp
  shuffler.hpp
  ga_shuffler.hpp
  printit.hpp
  task_manager.hpp
  random.hpp
  index_hash.hpp
  global_store.hpp
  global_vector.hpp
  DESTINATION include/gridpack/parallel
)

# -------------------------------------------------------------
# gridpack parallel library installation
# -------------------------------------------------------------
install(TARGETS 
  gridpack_parallel
  DESTINATION lib
)


# -------------------------------------------------------------
# TEST: greetings
# A simple MPI program (that does not use boost::test)
# -------------------------------------------------------------
add_executable(greetings test/greetings.cpp)
target_link_libraries(greetings ${target_libraries})

if (NOT USE_PROGRESS_RANKS)
  add_test(greetings_serial greetings)	
  set_tests_properties(greetings_serial
    PROPERTIES 
    PASS_REGULAR_EXPRESSION "I am process 0 of [1-9].*$"
    )
  set_tests_ldpath(greetings_serial)
endif()

add_test(greetings_parallel 
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} greetings ${MPIEXEC_POSTFLAGS})
set_tests_properties(greetings_parallel
  PROPERTIES 
  PASS_REGULAR_EXPRESSION "I am process 0 of [1-9].*$"
)
set_tests_ldpath(greetings_parallel)

# -------------------------------------------------------------
# TEST: task_test
# A simple program to test the task manager
# -------------------------------------------------------------
add_executable(task_test test/task_test.cpp ../timer/local_timer.cpp)
target_link_libraries(task_test ${target_libraries})

gridpack_add_run_test(task_test task_test "")

# -------------------------------------------------------------
# TEST: mpi_test
# A simple MPI test using boost::test
# -------------------------------------------------------------
add_executable(mpi_test test/mpi_test.cpp)
target_link_libraries(mpi_test ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})
add_dependencies(mpi_test external_build)
gridpack_add_unit_test(mpi_test mpi_test)

# -------------------------------------------------------------
# TEST: shuffle_test
# -------------------------------------------------------------
add_executable(shuffle_test test/shuffle_test.cpp)
target_link_libraries(shuffle_test ${target_libraries}
)
gridpack_add_unit_test(shuffle shuffle_test)

# -------------------------------------------------------------
# TEST: hash_test
# A simple program to test the task manager
# -------------------------------------------------------------
add_executable(hash_test test/hash_test.cpp)
target_link_libraries(hash_test ${target_libraries})

gridpack_add_unit_test(hash_test hash_test)

# -------------------------------------------------------------
# TEST: random_test
# A simple program to test the random number generator
# -------------------------------------------------------------
add_executable(random_test test/random_test.cpp)
target_link_libraries(random_test ${target_libraries})

gridpack_add_run_test(random_test random_test "")

# -------------------------------------------------------------
# TEST: store_test
# A simple program to test the global store module
# -------------------------------------------------------------
add_executable(store_test test/store_test.cpp)
target_link_libraries(store_test ${target_libraries})

gridpack_add_run_test(store_test store_test "")

# -------------------------------------------------------------
# TEST: vector_test
# A simple program to test the global vector module
# -------------------------------------------------------------
add_executable(vector_test test/vector_test.cpp)
target_link_libraries(vector_test ${target_libraries})

gridpack_add_run_test(vector_test vector_test "")
