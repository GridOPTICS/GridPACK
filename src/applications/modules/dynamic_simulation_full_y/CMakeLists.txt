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
# Last Change: 2015-04-27 10:30:17 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_components
    gridpack_partition
    gridpack_parallel
    gridpack_math
    gridpack_configuration
    gridpack_timer
    gridpack_powerflow_module
    gridpack_dsmatrix_components
    gridpack_ymatrix_components
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${GOSS_LIBRARY}
    ${APR_LIBRARY})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/base_classes
      ${CMAKE_CURRENT_SOURCE_DIR}/model_classes)

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_library(gridpack_dynamic_simulation_full_y_module
   dsf_app_module.cpp
   dsf_factory_module.cpp
   dsf_components_module.cpp
   dsf_generator_factory_module.cpp
   base_classes/base_generator_model.cpp
   base_classes/base_exciter_model.cpp
   base_classes/base_governor_model.cpp
   model_classes/classical.cpp
)

# -------------------------------------------------------------
# component serialization tests
# -------------------------------------------------------------
target_link_libraries(gridpack_dynamic_simulation_full_y_module
                      ${target_libraries})
   
# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  dsf_app_module.hpp
  dsf_factory_module.hpp
  dsf_components_module.hpp
  dsf_generator_factory_module.hpp
  base_classes/base_generator_model.hpp
  base_classes/base_exciter_model.hpp
  base_classes/base_governor_model.hpp
  model_classes/classical.hpp
  DESTINATION include/gridpack/applications/modules/dynamic_simulation_full_y
)

install(TARGETS 
  gridpack_dynamic_simulation_full_y_module
  DESTINATION lib
)

# -------------------------------------------------------------
# Create simple test that runs state estimation code
# -------------------------------------------------------------
add_executable(dsf_test test/dsf_test.cpp)
target_link_libraries(dsf_test gridpack_dynamic_simulation_full_y_module
                     ${target_libraries})

add_custom_target(dsf_test_input
 
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/test/input.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/test/IEEE_145bus_v23_PSLF.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/test/IEEE_145b_classical_model.dyr
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS 
  ${CMAKE_CURRENT_SOURCE_DIR}/test/input.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/test/IEEE_145bus_v23_PSLF.raw
  ${CMAKE_CURRENT_SOURCE_DIR}/test/IEEE_145b_classical_model.dyr
)

add_dependencies (dsf_test dsf_test_input)

gridpack_add_run_test(dsf_test dsf_test "input.xml")

