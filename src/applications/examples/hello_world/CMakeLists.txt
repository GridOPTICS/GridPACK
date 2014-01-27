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
# Last Change: 2013-08-09 11:41:42 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_components
    gridpack_factory
    gridpack_partition
    gridpack_parallel
    gridpack_math
    gridpack_timer
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${PETSC_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_executable(hello_world.x
   hw_app.cpp
   hw_components.cpp
   hw_main.cpp
)

target_link_libraries(hello_world.x ${target_libraries})

add_custom_target(hello_world.x.input

   COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_SOURCE_DIR}/10x10.raw
   ${CMAKE_CURRENT_BINARY_DIR} 

   DEPENDS
   ${CMAKE_CURRENT_SOURCE_DIR}/10x10.raw

)
add_dependencies(hello_world.x hello_world.x.input)
