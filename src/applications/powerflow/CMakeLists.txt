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
    gridpack_stream
    gridpack_partition
    gridpack_environment
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

if (GOSS_DIR)
  set(target_libraries
      ${target_libraries}
      gridpack_goss
      ${GOSS_LIBRARY}
      ${APR_LIBRARY})
endif()

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
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_14.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/powerflow/input_14.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_14.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/powerflow/input_14.xml"
  )

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_118.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/powerflow/input_118.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_118.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/powerflow/input_118.xml"
  )

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_polish.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/powerflow/input_polish.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_polish.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/powerflow/input_polish.xml"
  )

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_european.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/powerflow/input_european.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_european.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/powerflow/input_european.xml"
  )


add_custom_target(pf.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/IEEE14.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/gridpack.petscrc
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/IEEE118.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/Polish_model_v23.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/EuropeanOpenModel_v23.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_BINARY_DIR}/input_14.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE14.raw
  ${CMAKE_CURRENT_SOURCE_DIR}/gridpack.petscrc
  ${CMAKE_CURRENT_BINARY_DIR}/input_118.xml
  ${GRIDPACK_DATA_DIR}/raw/IEEE118.raw
  ${CMAKE_CURRENT_BINARY_DIR}/input_polish.xml
  ${GRIDPACK_DATA_DIR}/raw/Polish_model_v23.raw
  ${CMAKE_CURRENT_BINARY_DIR}/input_european.xml
  ${GRIDPACK_DATA_DIR}/raw/EuropeanOpenModel_v23.raw

)
add_dependencies(pf.x pf.x.input)

# -------------------------------------------------------------
# install as a sample application
# -------------------------------------------------------------

install(TARGETS pf.x DESTINATION bin)
  
# -------------------------------------------------------------
# Create simple test that runs powerflow code
# -------------------------------------------------------------
gridpack_add_run_test("powerflow" pf.x "input_14.xml")

