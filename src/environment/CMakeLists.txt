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
# Last Change: 2020-01-07 08:51:33 d3g096
# -------------------------------------------------------------

add_library(gridpack_environment 
  environment.cpp
)

target_link_libraries(gridpack_environment gridpack_math)
  
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# -------------------------------------------------------------
# gridpack header installation
# -------------------------------------------------------------
install(FILES
  environment.hpp
  DESTINATION include/gridpack/environment
)

# -------------------------------------------------------------
# gridpack environment library installation
# -------------------------------------------------------------
install(TARGETS 
  gridpack_environment
  DESTINATION lib
)
