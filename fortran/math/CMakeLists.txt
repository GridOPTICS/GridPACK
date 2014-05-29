# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May 12, 2014 by William A. Perkins
# Last Change: 2014-05-29 13:50:13 d3g096
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

add_executable(linear_solver_test
  linear_solver_test.F90
  $<TARGET_OBJECTS:gridpackf_math>
  $<TARGET_OBJECTS:gridpackf_parallel>
)
target_link_libraries(linear_solver_test
  ${GRIDPACK_LIBS}
)

