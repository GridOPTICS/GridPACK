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
# Last Change: 2014-02-05 08:10:28 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_components
    gridpack_factory
    gridpack_partition
    gridpack_parallel
    gridpack_math
    gridpack_configuration
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

add_executable(resistor_grid.x
   rg_app.cpp
   rg_components.cpp
   rg_main.cpp
)

target_link_libraries(resistor_grid.x ${target_libraries})

add_custom_target(resistor_grid.x.input

   COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_SOURCE_DIR}/small.raw
   ${CMAKE_CURRENT_BINARY_DIR} 

   COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
   ${CMAKE_CURRENT_BINARY_DIR} 

   DEPENDS
   ${CMAKE_CURRENT_SOURCE_DIR}/small.raw
   ${CMAKE_CURRENT_SOURCE_DIR}/input.xml

)
add_dependencies(resistor_grid.x resistor_grid.x.input)

# -------------------------------------------------------------
# install as a sample application
# -------------------------------------------------------------

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.install.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)

install(FILES 
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
  ${CMAKE_CURRENT_SOURCE_DIR}/rg_app.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/rg_components.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/rg_main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/rg_app.hpp  
  ${CMAKE_CURRENT_SOURCE_DIR}/rg_components.hpp  
  ${CMAKE_CURRENT_SOURCE_DIR}/rg_factory.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/rg_grid.F90
  ${CMAKE_CURRENT_BINARY_DIR}/input.xml
  ${CMAKE_CURRENT_BINARY_DIR}/small.raw
  DESTINATION doc/gridpack/example/resistor_grid
)
