# -------------------------------------------------------------
# file: GridPACK.cmake
#
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created June 10, 2013 by William A. Perkins
# Last Change: 2013-06-10 10:44:41 d3g096
# -------------------------------------------------------------

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
    FAIL_REGULAR_EXPRESSION "failure"
  )
endfunction(gridpack_add_serial_unit_test)

# -------------------------------------------------------------
# gridpack_add_parallel_unit_test
# -------------------------------------------------------------
function(gridpack_add_parallel_unit_test test_name test_program)
  set(the_test_name "${test_name}_parallel")
  add_test("${the_test_name}"
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${test_program} ${MPIEXEC_POSTFLAGS})
  set_tests_properties("${the_test_name}"
    PROPERTIES 
    PASS_REGULAR_EXPRESSION "No errors detected"
    FAIL_REGULAR_EXPRESSION "failure"
  )
endfunction(gridpack_add_parallel_unit_test)


# -------------------------------------------------------------
# gridpack_add_unit_test
#
# A way to consistently add both a serial and parallel test of the
# same executable
# -------------------------------------------------------------

function(gridpack_add_unit_test test_name test_program)
  gridpack_add_serial_unit_test("${test_name}" "${test_program}")
  if (MPIEXEC) 
    gridpack_add_parallel_unit_test("${test_name}" "${test_program}")
  endif ()
endfunction(gridpack_add_unit_test)
