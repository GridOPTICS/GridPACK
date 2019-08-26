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
# Created June 11, 2013 by William A. Perkins
# Last Change: 2018-10-19 13:44:01 d3g096
# -------------------------------------------------------------

set(target_libraries 
  ${Boost_LIBRARIES} 
)

# -------------------------------------------------------------
# gridpack header installation
# -------------------------------------------------------------
install(FILES 
  complex.hpp
  exception.hpp
  uncopyable.hpp
  named.hpp
  doxygen_dummy.hpp
  indexed.hpp
  sanity_interface.hpp
  string_utils.hpp
  DESTINATION include/gridpack/utilities
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

# -------------------------------------------------------------
# TEST: test_network
# -------------------------------------------------------------
add_executable(test_string test/test_string.cpp)
target_link_libraries(test_string ${target_libraries})
add_dependencies(test_string external_build)