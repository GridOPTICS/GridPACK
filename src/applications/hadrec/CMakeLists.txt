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
    gridpack_hadrec_module
    gridpack_dynamic_simulation_full_y_module
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

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()
if (PARMETIS_FOUND)
  include_directories(BEFORE ${PARMETIS_INCLUDE_DIR})
endif()

add_executable(hadrec.x
   hadrec_main.cpp
)

target_link_libraries(hadrec.x ${target_libraries})

# Put files necessary to run hadrec.x in binary directory.
# gridpack.petscrc is temporary -- it will be incorporated into
# input.xml

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_tamu500_step005.xml"
  COMMAND ${CMAKE_COMMAND}
  -D INPUT:PATH="${GRIDPACK_DATA_DIR}/input/hadrec/input_tamu500_step005.xml"
  -D OUTPUT:PATH="${CMAKE_CURRENT_BINARY_DIR}/input_tamu500_step005.xml"
  -D PKG:STRING="${GRIDPACK_MATSOLVER_PKG}"
  -P "${PROJECT_SOURCE_DIR}/cmake-modules/set_lu_solver_pkg.cmake"
  DEPENDS "${GRIDPACK_DATA_DIR}/input/hadrec/input_tamu500_step005.xml"
  )



add_custom_target(hadrec.x.input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/raw/ACTIVSg500_rmsmallgen_pslfv23.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GRIDPACK_DATA_DIR}/dyr/tamu_500bus_detail.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_BINARY_DIR}/input_tamu500_step005.xml
  ${GRIDPACK_DATA_DIR}/raw/ACTIVSg500_rmsmallgen_pslfv23.raw
  ${GRIDPACK_DATA_DIR}/dyr/tamu_500bus_detail.dyr

)
add_dependencies(hadrec.x hadrec.x.input)

# -------------------------------------------------------------
# install as a sample application
# -------------------------------------------------------------

install(TARGETS hadrec.x DESTINATION bin)
  
# -------------------------------------------------------------
# Create simple test that runs hadrec code
# -------------------------------------------------------------
#gridpack_add_run_test("hadrec" hadrec.x "input_tamu500_step005.xml")
