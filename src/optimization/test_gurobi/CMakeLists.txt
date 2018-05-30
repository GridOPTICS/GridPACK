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
# Last Change: 2014-01-28 13:13:49 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_components
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

add_executable(uc.x
   uc_app.cpp
   uc_components.cpp
   uc_main.cpp
)

target_link_libraries(uc.x ${target_libraries})

add_custom_target(uc.x.input

   COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_SOURCE_DIR}/uc_test.raw
   ${CMAKE_CURRENT_SOURCE_DIR}/gen.uc
   ${CMAKE_CURRENT_BINARY_DIR} 

   DEPENDS
   ${CMAKE_CURRENT_SOURCE_DIR}/uc_test.raw
   ${CMAKE_CURRENT_SOURCE_DIR}/gen.uc

)
add_dependencies(uc.x uc.x.input)


# -------------------------------------------------------------
# install as a sample application
# -------------------------------------------------------------

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.install.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)

#install(FILES 
#  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
#  ${CMAKE_CURRENT_SOURCE_DIR}/uc_app.cpp
#  ${CMAKE_CURRENT_SOURCE_DIR}/uc_app.hpp
#  ${CMAKE_CURRENT_SOURCE_DIR}/uc_components.cpp
#  ${CMAKE_CURRENT_SOURCE_DIR}/uc_components.hpp
#  ${CMAKE_CURRENT_SOURCE_DIR}/uc_factory.hpp
#  ${CMAKE_CURRENT_SOURCE_DIR}/uc_main.cpp
#  ${CMAKE_CURRENT_SOURCE_DIR}/UC_test.raw
#  DESTINATION doc/gridpack/example/uc
#)

