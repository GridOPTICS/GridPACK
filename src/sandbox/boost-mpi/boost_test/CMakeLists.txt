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
# Last Change: 2014-01-28 13:13:59 d3g096
# -------------------------------------------------------------

set(target_libraries
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(shuffle.x
   shuffle.cpp
)

add_executable(redistribute.x
   redistribute.cpp
)

target_link_libraries(shuffle.x ${target_libraries})
target_link_libraries(redistribute.x ${target_libraries})
