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

add_library(gridpackf_parser OBJECT
  parser_c.cpp
  parser_f.F90
)

#add_dependencies( gridpackf_factory $<TARGET_OBJECTS:gridpackf_parser> )

#add_executable(fortran_parser_test
#  factory_test.F90
#  $<TARGET_OBJECTS:gridpackf_parser>
#  $<TARGET_OBJECTS:gridpackf_network>
#  $<TARGET_OBJECTS:gridpackf_parallel>
#  $<TARGET_OBJECTS:gridpackf_component>
#)

#target_link_libraries(fortran_parser_test
#  ${MPI_Fortran_LIBRARIES}
#  ${GRIDPACK_LIBS}
#)

