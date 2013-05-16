# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May  7, 2013 by William A. Perkins
# Last Change: 2013-05-16 09:02:20 d3g096
# -------------------------------------------------------------

set(gridpack_math_sources)
set(target_libraries gridpack_parallel ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})

# implementation independent math sources
list(APPEND gridpack_math_sources
  vector.cpp
  vector_implementation.cpp
  matrix.cpp
  matrix_implementation.cpp
  implementation_visitor.cpp
)

# implementation specific sources and flags
if (PETSC_FOUND)

    include_directories (AFTER ${PETSC_INCLUDES})
    add_definitions (${PETSC_DEFINITIONS})

    list(APPEND gridpack_math_sources 
      petsc/petsc_exception.cpp
      petsc/math.cpp
      petsc/vector.cpp
      petsc/matrix.cpp
      petsc/petsc_vector_implementation.cpp
      petsc/petsc_matrix_implementation.cpp
    )  

    set(target_libraries ${PETSC_LIBRARIES} ${target_libraries})

else()

    message(FATAL_ERROR "A math implementation (e.g. PETSc) must be specified")

endif()

# -------------------------------------------------------------
# gridpack math library
# -------------------------------------------------------------
add_library(gridpack_math
  ${gridpack_math_sources}
)

# -------------------------------------------------------------
# gridpack math header installation
# -------------------------------------------------------------


# -------------------------------------------------------------
# vector test suite
# -------------------------------------------------------------
add_executable(vector_test test/vector_test.cpp)
target_link_libraries(vector_test gridpack_math ${target_libraries})
add_test(NAME vector_test_serial COMMAND $<TARGET_FILE:vector_test>)
add_test(NAME vector_test_parallel COMMAND 
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:vector_test> ${MPIEXEC_POSTFLAGS})
set_tests_properties(vector_test_serial vector_test_parallel
  PROPERTIES 
  PASS_REGULAR_EXPRESSION "No errors detected"
)

add_executable(matrix_test test/matrix_test.cpp)
target_link_libraries(matrix_test gridpack_math ${target_libraries})
add_test(NAME matrix_test_serial COMMAND $<TARGET_FILE:matrix_test>)
add_test(NAME matrix_test_parallel COMMAND 
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:matrix_test> ${MPIEXEC_POSTFLAGS})
set_tests_properties(matrix_test_serial matrix_test_parallel
  PROPERTIES 
  PASS_REGULAR_EXPRESSION "No errors detected"
)
