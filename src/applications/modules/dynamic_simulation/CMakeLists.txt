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
# Created January 19, 2015 by William A. Perkins, Bruce Palmer
# -------------------------------------------------------------

set(target_libraries
    gridpack_dsmatrix_components
    gridpack_ymatrix_components
    gridpack_powerflow_module
    gridpack_components
    gridpack_stream
    gridpack_partition
    gridpack_parallel
    gridpack_environment
    gridpack_math
    gridpack_configuration
    gridpack_timer
    ${PETSC_LIBRARIES}
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

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

add_library(gridpack_dynamic_simulation_module
  ds_app_module.cpp
  ds_factory_module.cpp
)

gridpack_set_library_version(gridpack_dynamic_simulation_module)

target_link_libraries(gridpack_dynamic_simulation_module
  gridpack_dsmatrix_components
  gridpack_powerflow_module
  )

# -------------------------------------------------------------
# component serialization tests
# -------------------------------------------------------------
#target_link_libraries(gridpack_dynamic_simulation_module ${target_libraries})

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  ds_app_module.hpp
  ds_factory_module.hpp
  DESTINATION include/gridpack/applications/modules/dynamic_simulation
)

install(TARGETS 
  gridpack_dynamic_simulation_module
  DESTINATION lib
)

# -------------------------------------------------------------
# Create simple test that runs state estimation code
# -------------------------------------------------------------
#add_executable(ds_test test/ds_test.cpp)
#target_link_libraries (ds_test gridpack_dynamic_simulation_module
#${target_libraries})

#add_custom_target(ds_test_input
#
#  COMMAND ${CMAKE_COMMAND} -E copy 
#  ${CMAKE_CURRENT_SOURCE_DIR}/test/input.xml
#  ${CMAKE_CURRENT_BINARY_DIR}

#  COMMAND ${CMAKE_COMMAND} -E copy 
#  ${CMAKE_CURRENT_SOURCE_DIR}/test/IEEE_145bus_v23_PSLF.raw
#  ${CMAKE_CURRENT_BINARY_DIR}

#  COMMAND ${CMAKE_COMMAND} -E copy 
#  ${CMAKE_CURRENT_SOURCE_DIR}/test/IEEE_145b_classical_model.dyr
#  ${CMAKE_CURRENT_BINARY_DIR}

#  DEPENDS
#  ${CMAKE_CURRENT_SOURCE_DIR}/test/input.xml
#  ${CMAKE_CURRENT_SOURCE_DIR}/test/IEEE_145bus_v23_PSLF.raw
#  ${CMAKE_CURRENT_SOURCE_DIR}/test/IEEE_145b_classical_model.dyr
#)

#add_dependencies (ds_test ds_test_input)

#gridpack_add_run_test(ds_test ds_test "input.xml")
