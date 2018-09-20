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
# Last Change: 2018-09-20 11:25:03 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_ymatrix_components
    gridpack_components
    gridpack_parallel
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(gridpack_pfmatrix_components
  pf_components.cpp
)

# -------------------------------------------------------------
# component serialization tests
# -------------------------------------------------------------
if (GRIDPACK_LIB_LINK_LIBRARIES)
  target_link_libraries(gridpack_pfmatrix_components
    gridpack_ymatrix_components
    )
endif()

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  pf_components.hpp  
  DESTINATION include/gridpack/applications/components/pf_matrix
)

install(TARGETS 
  gridpack_pfmatrix_components
  DESTINATION lib
)


