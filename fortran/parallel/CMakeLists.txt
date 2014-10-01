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
# Created May 12, 2014 by William A. Perkins
# Last Change: 2014-08-06 13:54:54 d3g096
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

add_executable(mpitest mpitest.F90)
target_link_libraries(mpitest ${MPI_Fortran_LIBRARIES})
