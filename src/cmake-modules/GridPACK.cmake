# -*- mode: cmake -*-
#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be found
#     in the LICENSE file in the top level directory of this distribution.
#
# -------------------------------------------------------------
# file: GridPACK.cmake
#
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created June 10, 2013 by William A. Perkins
# Last Change: 2016-07-15 10:20:33 d3g096
# -------------------------------------------------------------

# This is used to specify a time out for GridPACK unit tests. It's 5
# seconds by default, but may need to be longer on some platforms.
if (NOT GRIDPACK_TEST_TIMEOUT) 
  set (GRIDPACK_TEST_TIMEOUT 5 
    CACHE STRING "Time out for GridPACK unit tests.")
endif ()


# -------------------------------------------------------------
# gridpack_add_serial_unit_test
#
# This provides a way to consistly add a test that uses Boost::Test
# and can be executed on one processor without ${MPI_EXEC}.
# -------------------------------------------------------------
function(gridpack_add_serial_unit_test test_name test_program)
  set(the_test_name "${test_name}_serial")
  add_test("${the_test_name}" "${test_program}")
  set_tests_properties("${the_test_name}"
    PROPERTIES 
    PASS_REGULAR_EXPRESSION "No errors detected"
    FAIL_REGULAR_EXPRESSION "failure detected"
  )
endfunction(gridpack_add_serial_unit_test)

# -------------------------------------------------------------
# gridpack_add_serial_run_test
#
# This provides a way to consistly add a test that just runs a program
# on one processor without using ${MPI_EXEC}. Success or failure is
# based on the exit code.
# -------------------------------------------------------------
function(gridpack_add_serial_run_test test_name test_program test_input)
  set(the_test_name "${test_name}_serial")
  add_test("${the_test_name}" "${test_program}" ${test_input})
  set_tests_properties("${the_test_name}"
    PROPERTIES 
    TIMEOUT ${GRIDPACK_TEST_TIMEOUT}
  )
endfunction(gridpack_add_serial_run_test)

# -------------------------------------------------------------
# gridpack_add_parallel_unit_test
# -------------------------------------------------------------
function(gridpack_add_parallel_unit_test test_name test_target)
  set(the_test_name "${test_name}_parallel")
  if (TARGET ${test_target})
    get_property(test_program TARGET ${test_target} PROPERTY LOCATION)
    add_test("${the_test_name}"
      ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ${test_program} ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${the_test_name}"
      PROPERTIES 
      PASS_REGULAR_EXPRESSION "No errors detected"
      FAIL_REGULAR_EXPRESSION "failure detected"
      TIMEOUT ${GRIDPACK_TEST_TIMEOUT}
      )
  else() 
    message(FATAL_ERROR "gridpack_add_parallel_unit_test: target arguement not target")
  endif()
endfunction(gridpack_add_parallel_unit_test)

# -------------------------------------------------------------
# gridpack_add_parallel_run_test
#
# This provides a way to consistly add a test that just runs a program
# on multiple processors using ${MPI_EXEC}. Success or failure is
# based on the exit code.
# -------------------------------------------------------------
function(gridpack_add_parallel_run_test test_name test_target test_input)
  set(the_test_name "${test_name}_parallel")
  if (TARGET ${test_target})
    get_property(test_program TARGET ${test_target} PROPERTY LOCATION)
    add_test("${the_test_name}"
      ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ${test_program} ${MPIEXEC_POSTFLAGS} ${test_input})
    set_tests_properties("${the_test_name}"
      PROPERTIES 
      TIMEOUT ${GRIDPACK_TEST_TIMEOUT}
    )
  else()
    message(FATAL_ERROR "gridpack_add_parallel_run_test: target arguement not target")
  endif()
endfunction(gridpack_add_parallel_run_test)


# -------------------------------------------------------------
# gridpack_add_unit_test
#
# A way to consistently add both a serial and parallel test of the
# same executable
# -------------------------------------------------------------

function(gridpack_add_unit_test test_name test_program)
  if (NOT USE_PROGRESS_RANKS)
    gridpack_add_serial_unit_test("${test_name}" ${test_program})
  endif()
  if (MPIEXEC) 
    gridpack_add_parallel_unit_test("${test_name}" ${test_program})
  endif ()
endfunction(gridpack_add_unit_test)


# -------------------------------------------------------------
# gridpack_add_run_test
#
# A way to consistently add both a serial and parallel run test of the
# same executable
# -------------------------------------------------------------

function(gridpack_add_run_test test_name test_target test_input)
  if (NOT USE_PROGRESS_RANKS)
    gridpack_add_serial_run_test("${test_name}" ${test_target}  "${test_input}")
  endif()
  if (MPIEXEC) 
    gridpack_add_parallel_run_test("${test_name}" ${test_target}  "${test_input}")
  endif ()
endfunction(gridpack_add_run_test)
