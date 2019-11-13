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
# Last Change: 2019-01-04 08:39:01 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_components
    gridpack_stream
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

add_executable(ds_implicit.x
   ds_app.cpp
   ds_components.cpp
   ds_factory.cpp
   ds_main.cpp
)

target_link_libraries(ds_implicit.x ${target_libraries})

add_custom_target(ds_implicit.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/case9.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/case9_classical.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/case9.raw
  ${CMAKE_CURRENT_SOURCE_DIR}/case9_classical.dyr
)

add_dependencies(ds_implicit.x ds_implicit.x.input)

# -------------------------------------------------------------
# install as an example
# -------------------------------------------------------------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.install.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)

# install(FILES 
#   ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
#   ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
#   ${CMAKE_CURRENT_SOURCE_DIR}/case9.raw
#   ${CMAKE_CURRENT_SOURCE_DIR}/case9_classical.dyr
#   ${CMAKE_CURRENT_SOURCE_DIR}/ds_app.cpp
#   ${CMAKE_CURRENT_SOURCE_DIR}/ds_app.hpp
#   ${CMAKE_CURRENT_SOURCE_DIR}/ds_components.cpp
#   ${CMAKE_CURRENT_SOURCE_DIR}/ds_components.hpp
#   ${CMAKE_CURRENT_SOURCE_DIR}/ds_factory.cpp
#   ${CMAKE_CURRENT_SOURCE_DIR}/ds_factory.hpp
#   ${CMAKE_CURRENT_SOURCE_DIR}/ds_main.cpp
#   DESTINATION share/gridpack/example/dynamic_simulation_implicit
# )

# install(TARGETS ds_implicit.x DESTINATION bin)
  

