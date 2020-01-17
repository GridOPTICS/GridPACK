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
# Last Change: 2018-09-20 11:18:01 d3g096
# -------------------------------------------------------------

add_library(gridpack_stream
  input_stream.cpp
)

gridpack_set_library_version(gridpack_stream)
add_dependencies(gridpack_stream external_build)

if (GRIDPACK_LIB_LINK_LIBRARIES)
  target_link_libraries(gridpack_stream)
endif()

set(target_libraries
    gridpack_stream
    ${GA_LIBRARIES}
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(AFTER ${GA_INCLUDE_DIRS})

# -------------------------------------------------------------
# TEST: test_stream
# -------------------------------------------------------------

#add_executable(test_stream test/test_stream.cpp)
#target_link_libraries(test_stream ${target_libraries})

#gridpack_add_unit_test(stream_test test_stream)

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  input_stream.hpp
  DESTINATION include/gridpack/stream
)

install(TARGETS 
  gridpack_stream
  DESTINATION lib
)
