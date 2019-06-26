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
# Last Change: 2018-09-20 11:17:39 d3g096
# -------------------------------------------------------------

add_library(gridpack_appenvironment 
  app_environment.cpp
)

# For some reason this is necessary for Mac OS X, but messes things up on Linux
if (GRIDPACK_LIB_LINK_LIBRARIES)
  target_link_libraries(gridpack_appenvironment gridpack_math
    ${Boost_LIBRARIES} 
    ${GA_LIBRARIES} 
    ${MPI_CXX_LIBRARIES}
    )
endif()
  
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# -------------------------------------------------------------
# gridpack header installation
# -------------------------------------------------------------
install(FILES
  app_environment.hpp
  DESTINATION include/gridpack/app_environment
)

# -------------------------------------------------------------
# gridpack app_environment library installation
# -------------------------------------------------------------
install(TARGETS 
  gridpack_appenvironment
  DESTINATION lib
)
