# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May 12, 2014 by William A. Perkins
# Last Change: 2014-08-22 08:15:29 d3g096
# -------------------------------------------------------------

add_library(gridpackf_math OBJECT
  mathc.cpp
  mathf.F90
  vectorc.cpp
  vectorf.F90
  matrixc.cpp
  matrixf.F90
  linear_solverc.cpp
  linear_solverf.F90
)

add_dependencies(gridpackf_math 
  gridpackf_parallel
  gridpackf_configuration
)

add_executable(hello_math
  hello_math.F90
  $<TARGET_OBJECTS:gridpackf_math>
  $<TARGET_OBJECTS:gridpackf_parallel>
)
target_link_libraries(hello_math
  ${GRIDPACK_LIBS}
)

add_executable(vector_test
  vector_test.F90
  $<TARGET_OBJECTS:gridpackf_math>
  $<TARGET_OBJECTS:gridpackf_parallel>
)
target_link_libraries(vector_test
  ${GRIDPACK_LIBS}
)

add_executable(matrix_test
  matrix_test.F90
  $<TARGET_OBJECTS:gridpackf_math>
  $<TARGET_OBJECTS:gridpackf_parallel>
)
target_link_libraries(matrix_test
  ${GRIDPACK_LIBS}
)

add_custom_target(math_input

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/gridpack.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/gridpack.xml
)

add_executable(linear_solver_test
  linear_solver_test.F90
  $<TARGET_OBJECTS:gridpackf_math>
  $<TARGET_OBJECTS:gridpackf_configuration>
  $<TARGET_OBJECTS:gridpackf_parallel>
)
target_link_libraries(linear_solver_test
  ${GRIDPACK_LIBS}
)
add_dependencies(linear_solver_test math_input)

