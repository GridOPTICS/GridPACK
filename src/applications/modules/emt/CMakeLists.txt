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

add_library(gridpack_emt_module
  emt.cpp
  emtfactory.cpp
  emtnetwork.cpp
  emtutilfunctions.cpp
  cblock.cpp
  dblock.cpp
  base_classes/base_gen_model.cpp
  base_classes/base_exc_model.cpp
  base_classes/base_gov_model.cpp
  base_classes/base_load_model.cpp
  base_classes/base_plant_model.cpp
  base_classes/base_branch_model.cpp
  base_classes/fault.cpp	
  model_classes/gencls.cpp
  model_classes/genrou.cpp
  model_classes/constantimpedance.cpp
  model_classes/exdc1.cpp
  model_classes/ieeet1.cpp
  model_classes/wsieg1.cpp
  model_classes/tgov1.cpp
  model_classes/regca1.cpp
  model_classes/reeca1.cpp
  model_classes/repca1.cpp
  model_classes/gdform.cpp
  model_classes/lumpedline.cpp
  model_classes/transformer.cpp
  model_classes/sexs.cpp
)

target_link_libraries(gridpack_emt_module
  gridpack_ymatrix_components
  gridpack_partition
)

# -------------------------------------------------------------
# component serialization tests
# -------------------------------------------------------------
# target_link_libraries(gridpack_emt_module
#                       ${target_libraries})
   
# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  emt.hpp
  constants.hpp
  emtfactory.hpp
  emtnetwork.hpp
  emtutils.hpp
  cblock.hpp
  dblock.hpp
  emtutilfunctions.hpp
  DESTINATION include/gridpack/applications/modules/emt
)

install(FILES 
  base_classes/base_gen_model.hpp
  base_classes/base_exc_model.hpp
  base_classes/base_gov_model.hpp
  base_classes/base_load_model.hpp
  base_classes/base_plant_model.hpp
  base_classes/base_branch_model.hpp
  base_classes/fault.hpp
  DESTINATION include/gridpack/applications/modules/emt/base_classes
)

install(FILES
  model_classes/gencls.hpp
  model_classes/genrou.hpp
  model_classes/constantimpedance.hpp
  model_classes/exdc1.hpp
  model_classes/ieeet1.hpp
  model_classes/wsieg1.hpp
  model_classes/tgov1.hpp
  model_classes/regca1.hpp
  model_classes/reeca1.hpp
  model_classes/repca1.hpp
  model_classes/gdform.hpp
  model_classes/lumpedline.hpp
  model_classes/transformer.hpp
  model_classes/sexs.hpp
  DESTINATION include/gridpack/applications/modules/emt/model_classes
)

install(TARGETS 
  gridpack_emt_module
  DESTINATION lib
)
