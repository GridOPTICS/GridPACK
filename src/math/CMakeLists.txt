# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# Copyright (c) 2013 Battelle Memorial Institute
# Licensed under modified BSD License. A copy of this license can be
# found in the LICENSE file in the top level directory of this
# distribution.
# -------------------------------------------------------------
# Created May  7, 2013 by William A. Perkins
# Last Change: 2014-10-31 13:17:27 d3g096
# -------------------------------------------------------------
set(gridpack_math_sources)
set(target_libraries 
  gridpack_configuration 
  gridpack_parallel 
  gridpack_timer
  ${GA_LIBRARIES}
  ${Boost_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
)

# implementation independent math sources
list(APPEND gridpack_math_sources
  # matrix.cpp
  # matrix_implementation.cpp
  # linear_solver.cpp
  # linear_solver_implementation.cpp
  # linear_matrix_solver.cpp
  # linear_matrix_solver_implementation.cpp
  # basic_linear_matrix_solver_implementation.cpp
  # nonlinear_solver_interface.hpp
  # nonlinear_solver_interface.cpp
  # nonlinear_solver.cpp
  # nonlinear_solver_implementation.cpp
  # newton_raphson_solver.cpp
  # newton_raphson_solver_implementation.cpp 
  # dae_solver.cpp
  # dae_solver_implementation.cpp
  implementation_visitor.cpp
)

# implementation independent math header files
set(gridpack_math_headers
  basic_linear_matrix_solver_implementation.hpp
  dae_solver.hpp
  dae_solver_functions.hpp
  dae_solver_implementation.hpp
  implementation_visitor.hpp
  linear_matrix_solver.hpp
  linear_matrix_solver_implementation.hpp
  linear_solver.hpp
  linear_solver_implementation.hpp
  math.hpp
  matrix.hpp
  matrix_implementation.hpp
  newton_raphson_solver.hpp
  newton_raphson_solver_implementation.hpp
  nonlinear_solver.hpp
  nonlinear_solver_functions.hpp
  nonlinear_solver_implementation.hpp
  nonlinear_solver_interface.hpp
  vector.hpp
  vector_implementation.hpp
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
    include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# implementation specific sources and flags
if (PETSC_FOUND)

    include_directories (AFTER ${PETSC_INCLUDES})
    add_definitions (${PETSC_DEFINITIONS})

    if (PETSC_HAVE_SUPERLU_DIST) 
      add_definitions(-DPETSC_HAVE_SUPERLU)
    endif()

    list(APPEND gridpack_math_sources 
      petsc/petsc_exception.cpp
      petsc/petsc_configurable.cpp
      petsc/petsc_math.cpp
      petsc/petsc_vector.cpp
      petsc/petsc_vector_wrapper.cpp
      # petsc/petsc_matrix.cpp
      # petsc/petsc_matrix_implementation.cpp
      # petsc/petsc_matrix_operations.cpp
      # petsc/petsc_linear_solver.cpp
      # petsc/petsc_linear_solver_implementation.cpp
      # petsc/petsc_linear_matrix_solver.cpp
      # petsc/petsc_linear_matrix_solver_impl.cpp
      # petsc/petsc_nonlinear_solver_implementation.cpp
      # petsc/petsc_nonlinear_solver.cpp
      # petsc/petsc_dae_solver_implementation.cpp
      # petsc/petsc_dae_solver.cpp
    )  

    set(target_libraries ${PETSC_LIBRARIES} ${target_libraries})

    # PETSc-specific GridPACK configuration file(s) for test suites
    # gridpack.petsrc is temporary.  It's used to set the PETSc options in the
    # linear solver test.  It should not be necessary when a usable
    # configuration system is available.
    add_custom_target(math_test_input

      COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_SOURCE_DIR}/petsc/gridpack.petscrc
      ${CMAKE_CURRENT_BINARY_DIR}

      COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_SOURCE_DIR}/petsc/gridpack.xml
      ${CMAKE_CURRENT_BINARY_DIR}

      DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/petsc/gridpack.petscrc
      ${CMAKE_CURRENT_SOURCE_DIR}/petsc/gridpack.xml

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
# vector test suite
# -------------------------------------------------------------
add_executable(vector_test test/vector_test.cpp)
target_link_libraries(vector_test gridpack_math ${target_libraries})
gridpack_add_unit_test(vector vector_test)

# -------------------------------------------------------------
# matrix test suite (one each for dense and sparse)
# -------------------------------------------------------------
# add_executable(sparse_matrix_test test/matrix_test.cpp)
# target_link_libraries(sparse_matrix_test gridpack_math ${target_libraries})
# gridpack_add_unit_test(sparse_matrix sparse_matrix_test)

# add_executable(dense_matrix_test test/matrix_test.cpp)
# set_target_properties(dense_matrix_test
#   PROPERTIES
#   COMPILE_DEFINITIONS "TEST_DENSE=YES"
# )
# target_link_libraries(dense_matrix_test gridpack_math ${target_libraries})
# gridpack_add_unit_test(dense_matrix dense_matrix_test)

# -------------------------------------------------------------
# linear solver test suite
# -------------------------------------------------------------
# add_executable(linear_solver_test test/linear_solver_test.cpp)
# target_link_libraries(linear_solver_test gridpack_math ${target_libraries})
# gridpack_add_unit_test(linear_solver linear_solver_test)
# if (PETSC_FOUND)
#    add_dependencies(linear_solver_test math_test_input)
# endif()

# -------------------------------------------------------------
# nonlinear solver test suite
# -------------------------------------------------------------
# add_executable(nonlinear_solver_test test/nonlinear_solver_test.cpp)
# target_link_libraries(nonlinear_solver_test gridpack_math ${target_libraries})
# gridpack_add_unit_test(nonlinear_solver nonlinear_solver_test)
# if (PETSC_FOUND)
#    add_dependencies(nonlinear_solver_test math_test_input)
# endif()

# -------------------------------------------------------------
# DAE solver test suite
# -------------------------------------------------------------
# add_executable(dae_solver_test test/dae_solver_test.cpp)
# target_link_libraries(dae_solver_test gridpack_math ${target_libraries})
# gridpack_add_unit_test(dae_solver dae_solver_test)
# if (PETSC_FOUND)
#    add_dependencies(dae_solver_test math_test_input)
# endif()

  


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

