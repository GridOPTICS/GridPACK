# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created August 5, 2014 by Bruce Palmer
# Last Change: 2014-08-06 13:53:21 d3g096
# -------------------------------------------------------------

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(BEFORE ${GRIDPACK_INCLUDE_DIRS})

add_library(gridpackf_mapper OBJECT
  mapper_c.cpp
  mapper_f.F90
)

#add_dependencies( gridpackf_mapper $<TARGET_OBJECTS:gridpackf_network> )

add_executable(fortran_mapper_test
  mapper_test.F90
  mapper_component.F90
  mapper_factory.F90
  $<TARGET_OBJECTS:gridpackf_network>
  $<TARGET_OBJECTS:gridpackf_factory>
  $<TARGET_OBJECTS:gridpackf_mapper>
  $<TARGET_OBJECTS:gridpackf_component>
  $<TARGET_OBJECTS:gridpackf_parallel>
  $<TARGET_OBJECTS:gridpackf_math>
)

target_link_libraries(fortran_mapper_test
  ${MPI_Fortran_LIBRARIES}
  ${GRIDPACK_LIBS}
)
