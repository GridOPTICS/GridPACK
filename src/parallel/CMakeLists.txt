# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May  6, 2013 by William A. Perkins
# Last Change: 2013-05-07 13:57:25 d3g096
# -------------------------------------------------------------

add_library(gridpack_parallel 
  distributed.cpp
)

# -------------------------------------------------------------
# TEST: greetings
# A simple MPI program
# -------------------------------------------------------------
add_executable(greetings test/greetings.cpp)
target_link_libraries(greetings ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})

add_test(greetings_serial greetings)
add_test(greetings_parallel 
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} greetings ${MPIEXEC_POSTFLAGS})
set_tests_properties(greetings_serial greetings_parallel
  PROPERTIES 
  PASS_REGULAR_EXPRESSION "I am process 0 of [1-9].*$"
)

# -------------------------------------------------------------
# TEST: mpi_test
# A simple MPI test using boost::test
# -------------------------------------------------------------
add_executable(mpi_test test/mpi_test.cpp)
target_link_libraries(mpi_test ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})

add_test(mpi_test_serial mpi_test)
add_test(mpi_test_parallel 
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} mpi_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(mpi_test_serial mpi_test_parallel
  PROPERTIES 
  PASS_REGULAR_EXPRESSION "No errors detected"
)

