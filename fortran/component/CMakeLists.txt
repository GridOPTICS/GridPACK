#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be
#     found
#     in the LICENSE file in the top level directory of this distribution.
#
# -*- mode: cmake -*-
#
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May 12, 2014 by William A. Perkins
# Last Change: 2014-05-15 09:14:52 d3g096
# -------------------------------------------------------------

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(BEFORE ${GRIDPACK_INCLUDE_DIRS})

add_library(gridpackf_component OBJECT
  component_c.cpp
  component_f.F90
  data_collection_c.cpp
  data_collection_f.F90
  fortran_component.cpp
)

#add_executable(hello_parallel
#  hello.F90
#  $<TARGET_OBJECTS:gridpackf_component>
#)

#target_link_libraries(hello_parallel
#  ${GRIDPACK_LIBS}
#)

