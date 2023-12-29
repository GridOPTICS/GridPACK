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
# Last Change: 2023-12-13 07:00:00 d3g096
# -------------------------------------------------------------



# This is used to specify a time out for GridPACK unit tests. It's 5
# seconds by default, but may need to be longer on some platforms.
if (NOT GRIDPACK_TEST_TIMEOUT)
  set (GRIDPACK_TEST_TIMEOUT 120
    CACHE STRING "Time out for GridPACK unit tests.")
endif ()

# -------------------------------------------------------------
# set_tests_ldpath
# Use a macro to set each test's environment
# -------------------------------------------------------------
macro(set_tests_ldpath tname)
  if(GRIDPACK_TEST_LDPATH)
    if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
      set_tests_properties(${tname} PROPERTIES
        ENVIRONMENT DYLD_LIBRARY_PATH=${GRIDPACK_TEST_LDPATH}
        )
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
      set_tests_properties(${tname} PROPERTIES
        ENVIRONMENT LD_LIBRARY_PATH=${GRIDPACK_TEST_LDPATH}
        )
    endif()
  endif()
endmacro()

# -------------------------------------------------------------
# gridpack_set_library_version
# -------------------------------------------------------------
macro(gridpack_set_library_version lname)
  if(BUILD_SHARED_LIBS)
    set_target_properties(${lname}
      PROPERTIES
      SOVERSION "${GridPACK_VERSION_MAJOR}.${GridPACK_VERSION_MINOR}")
  endif()
endmacro()

# -------------------------------------------------------------
# gridpack_add_serial_unit_test
#
# This provides a way to consistly add a test that uses Boost::Test
# and can be executed on one processor without ${MPI_EXEC}.
# -------------------------------------------------------------
function(gridpack_add_serial_unit_test test_name test_target)
  set(the_test_name "${test_name}_serial")
  add_test("${the_test_name}" "${test_target}")
  set_tests_properties("${the_test_name}"
    PROPERTIES
    PASS_REGULAR_EXPRESSION "No errors detected"
    FAIL_REGULAR_EXPRESSION "failure detected"
    TIMEOUT ${GRIDPACK_TEST_TIMEOUT}
    )
  set_tests_ldpath("${the_test_name}")
endfunction(gridpack_add_serial_unit_test)

# -------------------------------------------------------------
# gridpack_add_serial_run_test
#
# This provides a way to consistly add a test that just runs a program
# on one processor without using ${MPI_EXEC}. Success or failure is
# based on the exit code.
# -------------------------------------------------------------
function(gridpack_add_serial_run_test test_name test_target test_input)
  set(the_test_name "${test_name}_serial")
  add_test("${the_test_name}" "${test_target}" ${test_input})
  set_tests_properties("${the_test_name}"
    PROPERTIES
    TIMEOUT ${GRIDPACK_TEST_TIMEOUT}
    )
  set_tests_ldpath("${the_test_name}")
endfunction(gridpack_add_serial_run_test)

# -------------------------------------------------------------
# gridpack_add_parallel_unit_test
# -------------------------------------------------------------
function(gridpack_add_parallel_unit_test test_name test_target)
  set(the_test_name "${test_name}_parallel")
  if (TARGET ${test_target})
    add_test("${the_test_name}"
      ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ${test_target} ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${the_test_name}"
      PROPERTIES
      PASS_REGULAR_EXPRESSION "No errors detected"
      FAIL_REGULAR_EXPRESSION "failure detected"
      TIMEOUT ${GRIDPACK_TEST_TIMEOUT}
      )
    set_tests_ldpath("${the_test_name}")
  else()
    message(FATAL_ERROR "gridpack_add_parallel_unit_test: target argument not target")
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
    add_test("${the_test_name}"
      ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ${test_target} ${MPIEXEC_POSTFLAGS} ${test_input})
    set_tests_properties("${the_test_name}"
      PROPERTIES
      TIMEOUT ${GRIDPACK_TEST_TIMEOUT}
      )
    set_tests_ldpath("${the_test_name}")
  else()
    message(FATAL_ERROR "gridpack_add_parallel_run_test: target argument not target")
  endif()
endfunction(gridpack_add_parallel_run_test)


# -------------------------------------------------------------
# gridpack_add_unit_test
#
# A way to consistently add both a serial and parallel test of the
# same executable
# -------------------------------------------------------------

function(gridpack_add_unit_test test_name test_target)
  if (NOT USE_PROGRESS_RANKS)
    gridpack_add_serial_unit_test("${test_name}" ${test_target})
  endif()
  if (MPIEXEC)
    gridpack_add_parallel_unit_test("${test_name}" ${test_target})
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


# -------------------------------------------------------------
# gridpack_set_lu_solver
#
# This takes a GridPACK input XML file, looks for any occurance of
# "superlu_dist" and changes it to "mumps" if SuperLU_DIST is not
# available in the build. The result is copied to the specified
# output.
# -------------------------------------------------------------
macro(gridpack_set_lu_solver inputxml outputxml)
add_custom_command(
  OUTPUT "${outputxml}"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${inputxml}"
  -D OUTPUT:PATH="${outputxml}"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${inputxml}"
)
endmacro()
