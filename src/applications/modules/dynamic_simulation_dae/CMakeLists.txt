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
# Last Change: 2022-09-16 07:49:18 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_ymatrix_components
    gridpack_components
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

set(ds_include_dirs
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/base_classes
    ${CMAKE_CURRENT_SOURCE_DIR}/model_classes)
include_directories(BEFORE ${ds_include_dirs})

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_library(gridpack_dynamic_simulation_dae_module
  dsim.cpp
  dsimfactory.cpp
  dsimnetwork.cpp
  base_classes/base_gen_model.cpp
  base_classes/base_exc_model.cpp
  base_classes/base_gov_model.cpp
  model_classes/classical_gen_model.cpp
  model_classes/genrou.cpp
  model_classes/gensal.cpp
  model_classes/esst1a.cpp
  model_classes/exdc1.cpp
  model_classes/wsieg1.cpp
  model_classes/GainBlockClass.cpp
  model_classes/BackLashClass.cpp
  model_classes/DBIntClass.cpp
)

target_link_libraries(gridpack_dynamic_simulation_dae_module
  gridpack_ymatrix_components
)

# -------------------------------------------------------------
# component serialization tests
# -------------------------------------------------------------
# target_link_libraries(gridpack_dynamic_simulation_dae_module
#                       ${target_libraries})
   
# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  dsim.hpp
  constants.hpp
  dsimfactory.hpp
  dsimnetwork.hpp
  dsimutils.hpp
  DESTINATION include/gridpack/applications/modules/dynamic_simulation_dae
)

install(FILES 
  base_classes/base_gen_model.hpp
  base_classes/base_exc_model.hpp
  base_classes/base_gov_model.hpp
  DESTINATION include/gridpack/applications/modules/dynamic_simulation_dae/base_classes
)

install(FILES 
  model_classes/classical_gen_model.hpp
  model_classes/genrou.hpp
  model_classes/gensal.hpp
  model_classes/exdc1.hpp
  model_classes/esst1a.hpp
  model_classes/wsieg1.hpp
  DESTINATION include/gridpack/applications/modules/dynamic_simulation_dae/model_classes
)

install(TARGETS 
  gridpack_dynamic_simulation_dae_module
  DESTINATION lib
)
