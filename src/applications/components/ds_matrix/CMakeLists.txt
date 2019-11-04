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
# Created May  6, 2013 by William A. Perkins
# Last Change: 2019-08-16 13:45:14 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_ymatrix_components
    gridpack_components
    gridpack_parallel
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(AFTER ${GA_INCLUDE_DIRS})

add_library(gridpack_dsmatrix_components
  ds_components.cpp
)

target_link_libraries(gridpack_dsmatrix_components
  gridpack_ymatrix_components
  )
gridpack_set_library_version(gridpack_dsmatrix_components)


# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  ds_components.hpp  
  DESTINATION include/gridpack/applications/components/ds_matrix
)

install(TARGETS 
  gridpack_dsmatrix_components
  DESTINATION lib
)


