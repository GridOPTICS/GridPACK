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
# Last Change: 2014-02-10 15:05:07 d3g096
# -------------------------------------------------------------

add_library(gridpack_parallel 
  environment.cpp
  distributed.cpp
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# -------------------------------------------------------------
# gridpack header installation
# -------------------------------------------------------------
install(FILES
  environment.hpp
  communicator.hpp
  parallel.hpp
  distributed.hpp
  shuffler.hpp
  printit.hpp
  task_manager.hpp
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
target_link_libraries(greetings gridpack_parallel 
  ${GA_LIBRARIES} ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})

add_test(greetings_serial greetings)
add_test(greetings_parallel 
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} greetings ${MPIEXEC_POSTFLAGS})
set_tests_properties(greetings_serial greetings_parallel
  PROPERTIES 
  PASS_REGULAR_EXPRESSION "I am process 0 of [1-9].*$"
)

# -------------------------------------------------------------
# TEST: task_test
# A simple program to test the task manager
# -------------------------------------------------------------
add_executable(task_test test/task_test.cpp)
target_link_libraries(task_test gridpack_parallel 
  ${GA_LIBRARIES} ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})

add_test(task_test_serial task_test)
add_test(task_test_parallel 
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} task_test ${MPIEXEC_POSTFLAGS})
#set_tests_properties(task_test_serial task_test_parallel
#  PROPERTIES 
#  PASS_REGULAR_EXPRESSION "I am process 0 of [1-9].*$"
#)

# -------------------------------------------------------------
# TEST: mpi_test
# A simple MPI test using boost::test
# -------------------------------------------------------------
add_executable(mpi_test test/mpi_test.cpp)
target_link_libraries(mpi_test ${GA_LIBRARIES} ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})
gridpack_add_unit_test(mpi_test mpi_test)


# -------------------------------------------------------------
# TEST: shuffle_test
# -------------------------------------------------------------
add_executable(shuffle_test test/shuffle_test.cpp)
target_link_libraries(shuffle_test ${GA_LIBRARIES} ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})
gridpack_add_unit_test(shuffle shuffle_test)
