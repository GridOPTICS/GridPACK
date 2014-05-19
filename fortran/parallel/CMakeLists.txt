# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May 12, 2014 by William A. Perkins
# Last Change: 2014-05-15 09:14:52 d3g096
# -------------------------------------------------------------

add_library(gridpackf_parallel OBJECT
  communicatorc.cpp
  communicatorf.F90
  parallelc.cpp
  parallelf.F90
)

add_executable(hello_parallel
  hello.F90
  $<TARGET_OBJECTS:gridpackf_parallel>
)
target_link_libraries(hello_parallel
  ${GRIDPACK_LIBS}
)

