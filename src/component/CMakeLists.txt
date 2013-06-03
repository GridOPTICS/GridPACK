# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May  6, 2013 by William A. Perkins
# Last Change: 2013-05-07 13:57:25 d3g096
# -------------------------------------------------------------

add_library(gridpack_components
  base_component.cpp
  data_collection.cpp
)

# -------------------------------------------------------------
# TEST: greetings
# A simple MPI program
# -------------------------------------------------------------
#!!add_executable(greetings test/greetings.cpp)
#!!target_link_libraries(greetings ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})

#!!add_test(NAME greetings_serial COMMAND $<TARGET_FILE:greetings>)
#!!add_test(NAME greetings_parallel COMMAND 
#!!  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:greetings> ${MPIEXEC_POSTFLAGS})
#!!set_tests_properties(greetings_serial greetings_parallel
#!!  PROPERTIES 
#!!  PASS_REGULAR_EXPRESSION "I am process 0 of [1-9].*$"
#!!)

# -------------------------------------------------------------
# TEST: mpi_test
# A simple MPI test using boost::test
# -------------------------------------------------------------
#!!add_executable(mpi_test test/mpi_test.cpp)
#!!target_link_libraries(mpi_test ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})

#!!add_test(NAME mpi_test_serial COMMAND $<TARGET_FILE:mpi_test>)
#!!add_test(NAME mpi_test_parallel COMMAND 
#!!  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mpi_test> ${MPIEXEC_POSTFLAGS})
#!!set_tests_properties(mpi_test_serial mpi_test_parallel
#!!  PROPERTIES 
#!!  PASS_REGULAR_EXPRESSION "No errors detected"
#!!)

