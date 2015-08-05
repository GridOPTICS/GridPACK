# -*- mode: cmake -*-
#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be
#     found
#     in the LICENSE file in the top level directory of this distribution.
#
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created June 11, 2013 by William A. Perkins
# Last Change: 2015-07-31 13:43:36 d3g096
# -------------------------------------------------------------

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(target_libraries 
  ${Boost_LIBRARIES} 
)

add_executable(variable_test variable_test.cpp variable.cpp)
target_link_libraries(variable_test ${target_libraries})

add_executable(expression_test expression_test.cpp variable.cpp)
target_link_libraries(expression_test ${target_libraries})


