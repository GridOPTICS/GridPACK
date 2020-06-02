# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# Copyright (c) 2013 Battelle Memorial Institute
# Licensed under modified BSD License. A copy of this license can be
# found in the LICENSE file in the top level directory of this
# distribution.
# -------------------------------------------------------------
# Created May  7, 2013 by William A. Perkins
# Last Change: 2020-06-02 06:48:10 d3g096
# -------------------------------------------------------------
set(gridpack_math_sources)
set(target_libraries 
  gridpack_configuration 
  gridpack_environment
  gridpack_math
  gridpack_parallel 
  gridpack_timer
  ${GA_LIBRARIES}
  ${Boost_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
  )

# implementation independent math sources
list(APPEND gridpack_math_sources
  vector_implementation.cpp
  implementation_visitor.cpp
  )

# implementation independent math header files
set(gridpack_math_headers
  basic_linear_matrix_solver_implementation.hpp
  dae_solver.hpp
  dae_solver_functions.hpp
  dae_solver_interface.hpp
  dae_solver_implementation.hpp
  complex_operators.hpp
  implementation_visitable.hpp
  implementation_visitor.hpp
  linear_matrix_solver.hpp
  linear_matrix_solver_implementation.hpp
  linear_matrix_solver_interface.hpp
  linear_solver.hpp
  linear_solver_implementation.hpp
  linear_solver_interface.hpp
  math.hpp
  matrix.hpp
  matrix_implementation.hpp
  matrix_interface.hpp
  matrix_storage_type.hpp
  newton_raphson_solver.hpp
  newton_raphson_solver_implementation.hpp
  nonlinear_solver_functions.hpp
  nonlinear_solver.hpp
  nonlinear_solver_implementation.hpp
  nonlinear_solver_interface.hpp
  vector.hpp
  vector_implementation.hpp
  vector_interface.hpp
  value_transfer.hpp
  numeric_type_check.hpp
  )

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# implementation specific sources and flags
if (PETSC_FOUND)

  include_directories (AFTER ${PETSC_INCLUDES})
  add_definitions (${PETSC_DEFINITIONS})

  list(APPEND gridpack_math_sources 
    petsc/petsc_exception.cpp
    petsc/petsc_configurable.cpp
    petsc/petsc_misc.cpp
    petsc/petsc_math.cpp
    petsc/petsc_vector.cpp
    petsc/petsc_vector_wrapper.cpp
    petsc/petsc_matrix.cpp
    petsc/ga_matrix.cpp
    petsc/petsc_matrix_operations.cpp
    petsc/petsc_matrix_wrapper.cpp
    petsc/petsc_matrix.cpp
    petsc/petsc_linear_solver.cpp
    petsc/petsc_linear_matrix_solver.cpp
    petsc/petsc_nonlinear_solver_implementation.cpp
    petsc/petsc_nonlinear_solver.cpp
    petsc/petsc_dae_solver.cpp
    )  

  set(target_libraries ${PETSC_LIBRARIES} ${target_libraries})

  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/gridpack.xml"
    COMMAND ${CMAKE_COMMAND}
    -D INPUT:PATH="${CMAKE_CURRENT_SOURCE_DIR}/petsc/gridpack.xml"
    -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/gridpack.xml"
    -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
    -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/petsc/gridpack.xml"
    )
  # PETSc-specific GridPACK configuration file(s) for test suites
  # gridpack.petsrc is temporary.  It's used to set the PETSc options in the
  # linear solver test.  It should not be necessary when a usable
  # configuration system is available.
  add_custom_target(math_test_input

    COMMAND ${CMAKE_COMMAND} -E copy 
    ${CMAKE_CURRENT_SOURCE_DIR}/petsc/gridpack.petscrc
    ${CMAKE_CURRENT_BINARY_DIR}

    COMMAND ${CMAKE_COMMAND} -E copy 
    ${CMAKE_CURRENT_SOURCE_DIR}/test/bouncing_ball_test.xml
    ${CMAKE_CURRENT_BINARY_DIR}

    DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/petsc/gridpack.petscrc
    ${CMAKE_CURRENT_SOURCE_DIR}/test/bouncing_ball_test.xml
    ${CMAKE_CURRENT_BINARY_DIR}/gridpack.xml
    )

else()

  message(FATAL_ERROR "A math implementation (e.g. PETSc) must be specified")

endif()

# -------------------------------------------------------------
# gridpack math library
# -------------------------------------------------------------
add_library(gridpack_math
  ${gridpack_math_sources}
  )

gridpack_set_library_version(gridpack_math)
add_dependencies(gridpack_math external_build)

target_link_libraries(gridpack_math 
  gridpack_configuration
  gridpack_timer
  gridpack_parallel
  ${PETSC_LIBRARIES}
  )

# -------------------------------------------------------------
# a small test (not a unit test)
# -------------------------------------------------------------
# add_executable(small_matrix_solve small_matrix_solve.cpp)
# target_link_libraries(small_matrix_solve  gridpack_math ${target_libraries})
# add_custom_target(small_matrix_solve_input

#   COMMAND ${CMAKE_COMMAND} -E copy 
#   ${CMAKE_CURRENT_SOURCE_DIR}/small_matrix_solve.xml
#   ${CMAKE_CURRENT_BINARY_DIR}

#   DEPENDS
#   ${CMAKE_CURRENT_SOURCE_DIR}/small_matrix_solve.xml
# )
# add_dependencies(small_matrix_solve small_matrix_solve_input)

# -------------------------------------------------------------
# matrix inverter program (not a unit test)
# -------------------------------------------------------------
add_executable(matrix_inverse matrix_inverse.cpp)
target_link_libraries(matrix_inverse  gridpack_math ${target_libraries})
add_custom_target(matrix_inverse_input

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/A.matrix
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/A.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/A.matrix
  ${CMAKE_CURRENT_SOURCE_DIR}/A.xml
  )
add_dependencies(matrix_inverse matrix_inverse_input)

# -------------------------------------------------------------
# numeric test suite
# -------------------------------------------------------------
add_executable(numeric_test test/numeric_test.cpp)
gridpack_add_serial_unit_test(numeric numeric_test)


# -------------------------------------------------------------
# vector test suite
# -------------------------------------------------------------
add_executable(complex_vector_test test/vector_test.cpp)
target_link_libraries(complex_vector_test gridpack_math ${target_libraries})
gridpack_add_unit_test(complex_vector complex_vector_test)

add_executable(real_vector_test test/vector_test.cpp)
set_target_properties(real_vector_test
  PROPERTIES
  COMPILE_DEFINITIONS "TEST_REAL=YES"
  )
target_link_libraries(real_vector_test gridpack_math ${target_libraries})
gridpack_add_unit_test(real_vector real_vector_test)


# -------------------------------------------------------------
# matrix test suite (one each for dense and sparse)
# -------------------------------------------------------------
add_executable(complex_sparse_matrix_test test/matrix_test.cpp)
target_link_libraries(complex_sparse_matrix_test gridpack_math ${target_libraries})
gridpack_add_unit_test(complex_sparse_matrix complex_sparse_matrix_test)

add_executable(real_sparse_matrix_test test/matrix_test.cpp)
set_target_properties(real_sparse_matrix_test
  PROPERTIES
  COMPILE_DEFINITIONS "TEST_REAL=YES"
  )
target_link_libraries(real_sparse_matrix_test gridpack_math ${target_libraries})
gridpack_add_unit_test(real_sparse_matrix real_sparse_matrix_test)

add_executable(complex_dense_matrix_test test/matrix_test.cpp)
set_target_properties(complex_dense_matrix_test
  PROPERTIES
  COMPILE_DEFINITIONS "TEST_DENSE=YES"
  )
target_link_libraries(complex_dense_matrix_test gridpack_math ${target_libraries})
gridpack_add_unit_test(complex_dense_matrix complex_dense_matrix_test)

add_executable(real_dense_matrix_test test/matrix_test.cpp)
set_target_properties(real_dense_matrix_test
  PROPERTIES
  COMPILE_DEFINITIONS "TEST_DENSE=YES;TEST_REAL=YES"
  )
target_link_libraries(real_dense_matrix_test gridpack_math ${target_libraries})
gridpack_add_unit_test(real_dense_matrix real_dense_matrix_test)

# -------------------------------------------------------------
# matrix transpose test
# -------------------------------------------------------------

add_executable(complex_sparse_transpose_test test/matrix_transpose_test.cpp)
target_link_libraries(complex_sparse_transpose_test gridpack_math ${target_libraries})
gridpack_add_unit_test(complex_sparse_matrix_transpose complex_sparse_transpose_test)

add_executable(real_sparse_transpose_test test/matrix_transpose_test.cpp)
set_target_properties(real_sparse_transpose_test
  PROPERTIES
  COMPILE_DEFINITIONS "TEST_REAL=YES"
  )
target_link_libraries(real_sparse_transpose_test gridpack_math ${target_libraries})
gridpack_add_unit_test(real_sparse_matrix_transpose real_sparse_transpose_test)

add_executable(complex_dense_transpose_test test/matrix_transpose_test.cpp)
set_target_properties(complex_dense_transpose_test
  PROPERTIES
  COMPILE_DEFINITIONS "TEST_DENSE=YES"
  )
target_link_libraries(complex_dense_transpose_test gridpack_math ${target_libraries})
gridpack_add_unit_test(complex_dense_matrix_transpose complex_dense_transpose_test)

add_executable(real_dense_transpose_test test/matrix_transpose_test.cpp)
set_target_properties(real_dense_transpose_test
  PROPERTIES
  COMPILE_DEFINITIONS "TEST_DENSE=YES;TEST_REAL=YES"
  )
target_link_libraries(real_dense_transpose_test gridpack_math ${target_libraries})
gridpack_add_unit_test(real_dense_matrix_transpose real_dense_transpose_test)


# -------------------------------------------------------------
# linear solver test suite
# -------------------------------------------------------------
add_executable(real_linear_solver_test test/linear_solver_test.cpp)
target_link_libraries(real_linear_solver_test gridpack_math ${target_libraries})
gridpack_add_unit_test(real_linear_solver real_linear_solver_test)
if (PETSC_FOUND)
  add_dependencies(real_linear_solver_test math_test_input)
endif()

add_executable(complex_linear_solver_test test/complex_linear_solver_test.cpp)
target_link_libraries(complex_linear_solver_test gridpack_math ${target_libraries})
gridpack_add_unit_test(complex_linear_solver complex_linear_solver_test)
if (PETSC_FOUND)
  add_dependencies(complex_linear_solver_test math_test_input)
endif()

# -------------------------------------------------------------
# nonlinear solver test suite
# -------------------------------------------------------------
add_executable(complex_nonlinear_solver_test test/nonlinear_solver_test.cpp)
target_link_libraries(complex_nonlinear_solver_test gridpack_math ${target_libraries})
gridpack_add_unit_test(complex_nonlinear_solver complex_nonlinear_solver_test)
if (PETSC_FOUND)
  add_dependencies(complex_nonlinear_solver_test math_test_input)
endif()

add_executable(real_nonlinear_solver_test test/nonlinear_solver_test.cpp)
set_target_properties(real_nonlinear_solver_test
  PROPERTIES
  COMPILE_DEFINITIONS "TEST_REAL=YES"
  )
target_link_libraries(real_nonlinear_solver_test gridpack_math ${target_libraries})
gridpack_add_unit_test(real_nonlinear_solver real_nonlinear_solver_test)
if (PETSC_FOUND)
  add_dependencies(real_nonlinear_solver_test math_test_input)
endif()


# -------------------------------------------------------------
# DAE solver test suite
# -------------------------------------------------------------
add_executable(complex_dae_solver_test test/dae_solver_test.cpp)
target_link_libraries(complex_dae_solver_test gridpack_math ${target_libraries})
gridpack_add_unit_test(complex_dae_solver complex_dae_solver_test)
if (PETSC_FOUND)
  add_dependencies(complex_dae_solver_test math_test_input)
endif()

add_executable(real_dae_solver_test test/dae_solver_test.cpp)
set_target_properties(real_dae_solver_test
  PROPERTIES
  COMPILE_DEFINITIONS "TEST_REAL=YES"
  )
target_link_libraries(real_dae_solver_test gridpack_math ${target_libraries})
gridpack_add_unit_test(real_dae_solver real_dae_solver_test)
if (PETSC_FOUND)
  add_dependencies(real_dae_solver_test math_test_input)
endif()

# A simple serial example utilizing DAE events
add_executable(bouncing_ball_test test/bouncing_ball_test.cpp)
target_link_libraries(bouncing_ball_test gridpack_math ${target_libraries})


# -------------------------------------------------------------
# Test GA-based Dense Matrix
# -------------------------------------------------------------
if (PETSC_FOUND)
  add_executable(petsc_ga_matrix_test test/petsc_ga_matrix.cpp) 
  target_link_libraries(petsc_ga_matrix_test gridpack_math ${target_libraries})
  gridpack_add_unit_test(petsc_ga_matrix petsc_ga_matrix_test)
endif()

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  ${gridpack_math_headers}
  DESTINATION include/gridpack/math
  )

install(TARGETS 
  gridpack_math
  DESTINATION lib
  )

