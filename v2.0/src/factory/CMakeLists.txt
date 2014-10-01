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
# Last Change: 2014-01-08 09:01:55 d3g096
# -------------------------------------------------------------

add_library(gridpack_factory
  base_factory.cpp
)

include_directories(AFTER ${GA_INCLUDE_DIRS})

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

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  base_factory.hpp
  DESTINATION include/gridpack/factory
)

install(TARGETS 
  gridpack_factory
  DESTINATION lib
)
