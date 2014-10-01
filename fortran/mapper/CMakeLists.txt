#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be
#     found
#     in the LICENSE file in the top level directory of this distribution.
#
# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created August 5, 2014 by Bruce Palmer
# Last Change: 2014-09-29 12:30:17 d3g096
# -------------------------------------------------------------

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(BEFORE ${GRIDPACK_INCLUDE_DIRS})

add_library(gridpackf_mapper OBJECT
  mapper_c.cpp
  full_map_f.F90
  bus_vector_map_f.F90
  gen_matrix_map_f.F90
  gen_vector_map_f.F90
)

add_dependencies(gridpackf_mapper 
  gridpackf_network
)

add_executable(fortran_mapper_test
  mapper_test.F90
  mapper_component.F90
  mapper_factory.F90
  $<TARGET_OBJECTS:gridpackf_network>
  $<TARGET_OBJECTS:gridpackf_factory>
  $<TARGET_OBJECTS:gridpackf_component>
  $<TARGET_OBJECTS:gridpackf_mapper>
  $<TARGET_OBJECTS:gridpackf_parallel>
  $<TARGET_OBJECTS:gridpackf_math>
)

target_link_libraries(fortran_mapper_test
  ${MPI_Fortran_LIBRARIES}
  ${GRIDPACK_LIBS}
)
