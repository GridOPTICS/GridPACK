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
# Last Change: 2014-08-06 13:24:48 d3g096
# -------------------------------------------------------------

set(target_libraries
    ${GA_LIBRARIES}
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

if (PARMETIS_FOUND)
  include_directories(BEFORE ${PARMETIS_INCLUDE_DIR})
  set(target_libraries 
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${target_libraries}
    )
endif()

set(target_libraries
    gridpack_partition
    gridpack_components
    gridpack_parallel
    gridpack_timer
    ${target_libraries})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()
# -------------------------------------------------------------
# TEST: parser_test
# A simple program to check input functionality
# -------------------------------------------------------------
#add_executable(PTI23_test test/PTI23_test.cpp)
#target_link_libraries(PTI23_test ${target_libraries})

add_executable(parser_test test/parser_test.cpp)
target_link_libraries(parser_test ${target_libraries})

add_custom_target(test_parser_input

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/test/test.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/test/parser_data.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GridPACK_SOURCE_DIR}/applications/powerflow/IEEE14.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/test/test.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/test/parser_data.raw
  ${GridPACK_SOURCE_DIR}/applications/powerflow/IEEE14.raw
)
add_dependencies(parser_test test_parser_input)
#add_dependencies(PTI23_test test_parser_input)

add_test(test_parser_serial parser_test)


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
  dictionary.hpp
  PTI23_parser.hpp
  hash_distr.hpp
  DESTINATION include/gridpack/parser
)
