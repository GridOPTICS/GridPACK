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
# Last Change: 2013-07-04 08:13:01 d3g096
# -------------------------------------------------------------

set(target_libraries
  gridpack_timer
  ${GA_LIBRARIES} 
  ${Boost_LIBRARIES} 
  ${MPI_CXX_LIBRARIES}
)

set(target_libraries 
  gridpack_components gridpack_parallel gridpack_math
  ${target_libraries} 
)

include_directories(AFTER ${GA_INCLUDE_DIRS})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# -------------------------------------------------------------
# TEST: greetings
# A simple Serial IO test program
# -------------------------------------------------------------
add_executable(test_serial_io test/test_serial_io.cpp)
target_link_libraries(test_serial_io ${target_libraries})
#
#!!add_test(test_network_serial_iotest_serial_io)
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

