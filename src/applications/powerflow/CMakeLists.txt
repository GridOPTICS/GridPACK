# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May  6, 2013 by William A. Perkins
# Last Change: 2013-10-07 15:12:17 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_components
    gridpack_factory
    gridpack_partition
    gridpack_parallel
    gridpack_math
    gridpack_configuration
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${PETSC_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_library(powerflow_app
  pf_components.cpp
  pf_factory.cpp
)

add_executable(powerflow.x
   pf_app.cpp
   pf_components.cpp
   pf_factory.cpp
   pf_main.cpp
)

target_link_libraries(powerflow.x ${target_libraries})

# Put files necessary to run powerflow.x in binary directory.
# gridpack.petscrc is temporary -- it will be incorporated into
# input.xml

add_custom_target(powerflow.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/IEEE14.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/gridpack.petscrc
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/118.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/118.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/IEEE14.raw
  ${CMAKE_CURRENT_SOURCE_DIR}/gridpack.petscrc
  ${CMAKE_CURRENT_SOURCE_DIR}/118.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/118.raw

)
add_dependencies(powerflow.x powerflow.x.input)

# -------------------------------------------------------------
# TEST: pf_network_test
# -------------------------------------------------------------
add_executable(pf_network_test test/pf_network_test.cpp)
target_link_libraries(pf_network_test powerflow_app ${target_libraries})
gridpack_add_unit_test(pf_network_test pf_network_test)

