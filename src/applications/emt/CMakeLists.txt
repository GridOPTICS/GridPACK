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
# Last Change: 2023-11-30 09:42:32 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_emt_module
    gridpack_powerflow_module
    gridpack_pfmatrix_components
    gridpack_ymatrix_components
    gridpack_components
    gridpack_partition
    gridpack_math
    gridpack_configuration
    gridpack_timer
    gridpack_parallel
    gridpack_block_parsers
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

include_directories(BEFORE
 ${CMAKE_CURRENT_SOURCE_DIR}/../modules/emt/model_classes)
include_directories(BEFORE
 ${CMAKE_CURRENT_SOURCE_DIR}/../modules/emt/base_classes)
include_directories(BEFORE
 ${CMAKE_CURRENT_SOURCE_DIR}/../modules/emt)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_executable(emt.x
   emt_app.cpp
)

target_link_libraries(emt.x ${target_libraries})

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/emt/input_emt.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/emt/input_emt.xml"
  )

add_custom_target(emt.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/case9.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENCLS.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENSAL.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU_ESST1A.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENSAL_ESST1A_WSIEG1.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU_EXDC1.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU_EXDC1_WSIEG1.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/9b3g.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/petscoptions/emt/.petscrc
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_BINARY_DIR}/input.xml
  ${GRIDPACK_DATA_DIR}/raw/case9.raw
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENCLS.dyr
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU.dyr
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU_ESST1A_WSIEG1.dyr
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENSAL.dyr
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU_ESST1A.dyr
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU_EXDC1.dyr
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU_EXDC1_WSIEG1.dyr
  ${GRIDPACK_DATA_DIR}/dyr/9b3g.dyr
)

add_dependencies(emt.x emt.x.input)

# -------------------------------------------------------------
# install as an example
# -------------------------------------------------------------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.install.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)

install(FILES 
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
  ${CMAKE_CURRENT_BINARY_DIR}/input.xml
  ${GRIDPACK_DATA_DIR}/raw/case9.raw
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENCLS.dyr
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU.dyr
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU_ESST1A_WSIEG1.dyr
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENSAL.dyr
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU_ESST1A.dyr
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU_EXDC1.dyr
  ${GRIDPACK_DATA_DIR}/dyr/case9_GENROU_EXDC1_WSIEG1.dyr
  ${GRIDPACK_DATA_DIR}/petscoptions/emt/.petscrc
  ${GRIDPACK_DATA_DIR}/dyr/9b3g.dyr
  ${CMAKE_CURRENT_SOURCE_DIR}/emt_app.cpp
  DESTINATION share/gridpack/example/emt
)

install(TARGETS emt.x DESTINATION bin)

# -------------------------------------------------------------
# run application as test
# -------------------------------------------------------------
gridpack_add_run_test("EMT" emt.x input.xml)

