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
# Last Change: 2017-12-07 08:54:26 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_powerflow_module
    gridpack_pfmatrix_components
    gridpack_ymatrix_components
    gridpack_components
    gridpack_partition
    gridpack_math
    gridpack_configuration
    gridpack_timer
    gridpack_parallel
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    )

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()
if (PARMETIS_FOUND)
  include_directories(BEFORE ${PARMETIS_INCLUDE_DIR})
endif()

add_executable(pf.x
   pf_main.cpp
)

target_link_libraries(pf.x ${target_libraries})

# Put files necessary to run pf.x in binary directory.
# gridpack.petscrc is temporary -- it will be incorporated into
# input.xml

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${CMAKE_CURRENT_SOURCE_DIR}/input.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/input.xml"
  )

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/118.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${CMAKE_CURRENT_SOURCE_DIR}/118.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/118.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/118.xml"
  )


add_custom_target(pf.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/IEEE14.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/gridpack.petscrc
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/118.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_BINARY_DIR}/input.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/IEEE14.raw
  ${CMAKE_CURRENT_SOURCE_DIR}/gridpack.petscrc
  ${CMAKE_CURRENT_BINARY_DIR}/118.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/118.raw

)
add_dependencies(pf.x pf.x.input)

# -------------------------------------------------------------
# install as a sample application
# -------------------------------------------------------------

# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.install.in
#   ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)

# install(FILES 
#   ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
#   ${CMAKE_CURRENT_SOURCE_DIR}/pf_main.cpp
#   ${CMAKE_CURRENT_BINARY_DIR}/input.xml
#   ${CMAKE_CURRENT_SOURCE_DIR}/IEEE14.raw
#   ${CMAKE_CURRENT_SOURCE_DIR}/gridpack.petscrc
#   ${CMAKE_CURRENT_BINARY_DIR}/118.xml
#   ${CMAKE_CURRENT_SOURCE_DIR}/118.raw
#   DESTINATION doc/gridpack/example/powerflow
# )

install(TARGETS pf.x DESTINATION bin)
  
# -------------------------------------------------------------
# Create simple test that runs powerflow code
# -------------------------------------------------------------
gridpack_add_run_test("powerflow" pf.x "input.xml")

