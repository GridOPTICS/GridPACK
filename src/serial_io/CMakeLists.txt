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
# Last Change: 2015-09-17 10:06:10 d3g096
# -------------------------------------------------------------

set(target_libraries
  gridpack_timer
  ${GA_LIBRARIES} 
  ${Boost_LIBRARIES} 
  ${MPI_CXX_LIBRARIES}
)

if (GOSS_DIR)
  set(target_libraries
      ${APR_LIBRARY}
      ${GOSS_LIBRARY}
      ${target_libraries})
endif()

set(target_libraries 
  gridpack_components gridpack_parallel gridpack_math
  gridpack_configuration
  ${target_libraries} 
)

include_directories(AFTER ${GA_INCLUDE_DIRS})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# -------------------------------------------------------------
# TEST: greetings
# A simple Serial IO test program
# -------------------------------------------------------------
add_executable(test_serial_io test/test_serial_io.cpp)
target_link_libraries(test_serial_io ${target_libraries})
if (GOSS_DIR)
  add_executable(test_goss test/test_goss.cpp)
  target_link_libraries(test_goss ${target_libraries} ${GOSS_LIBRARY} ${APR_LIBRARY})

  add_custom_target(test_goss_input
 
    COMMAND ${CMAKE_COMMAND} -E copy 
    ${CMAKE_CURRENT_SOURCE_DIR}/test/input_goss.xml
    ${CMAKE_CURRENT_BINARY_DIR}

    DEPENDS 
    ${CMAKE_CURRENT_SOURCE_DIR}/test/input_goss.xml
    )

  add_dependencies(test_goss test_goss_input)
endif()

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  serial_io.hpp
  DESTINATION include/gridpack/serial_io
)
