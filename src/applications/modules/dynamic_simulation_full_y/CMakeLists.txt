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
# Last Change: 2019-08-16 13:45:12 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_dsmatrix_components
    gridpack_ymatrix_components
    gridpack_powerflow_module
    gridpack_components
    gridpack_stream
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

if (GOSS_DIR)
  set(target_libraries
      ${target_libraries}
      gridpack_goss
      ${GOSS_LIBRARY}
      ${APR_LIBRARY})
endif()

if (HELICS_DIR)
  set(target_libraries
      ${target_libraries}
      ${HELICS_LIBRARY})
endif()

set(ds_include_dirs
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/base_classes
    ${CMAKE_CURRENT_SOURCE_DIR}/model_classes)
include_directories(BEFORE ${ds_include_dirs})

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_library(gridpack_dynamic_simulation_full_y_module
  dsf_app_module.cpp
  dsf_factory.cpp
  dsf_components.cpp
  generator_factory.cpp
  load_factory.cpp
  relay_factory.cpp
  base_classes/base_generator_model.cpp
  base_classes/base_exciter_model.cpp
  base_classes/base_pss_model.cpp
  base_classes/base_governor_model.cpp
  base_classes/base_relay_model.cpp
  base_classes/base_load_model.cpp
  model_classes/classical.cpp
  model_classes/gensal.cpp
  model_classes/exdc1.cpp
  model_classes/wsieg1.cpp
  model_classes/GainBlockClass.cpp
  model_classes/BackLashClass.cpp
  model_classes/DBIntClass.cpp
  model_classes/genrou.cpp
  model_classes/esst4b.cpp
  model_classes/esst1a.cpp
  model_classes/wshygp.cpp
  model_classes/ggov1.cpp
  model_classes/lvshbl.cpp
  model_classes/frqtpat.cpp
  model_classes/distr1.cpp
  model_classes/acmotor.cpp
  model_classes/ieel.cpp
  model_classes/motorw.cpp
  model_classes/psssim.cpp
)

gridpack_set_library_version(gridpack_dynamic_simulation_full_y_module)

target_link_libraries(gridpack_dynamic_simulation_full_y_module
  gridpack_dsmatrix_components
  gridpack_powerflow_module
  )

# -------------------------------------------------------------
# component serialization tests
# -------------------------------------------------------------
# target_link_libraries(gridpack_dynamic_simulation_full_y_module
#                       ${target_libraries})
   
# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  dsf_app_module.hpp
  dsf_components.hpp
  dsf_factory.hpp
  relay_factory.hpp
  generator_factory.hpp
  load_factory.hpp
  base_classes/base_generator_model.hpp
  base_classes/base_exciter_model.hpp
  base_classes/base_pss_model.hpp
  base_classes/base_governor_model.hpp
  base_classes/base_relay_model.hpp
  base_classes/base_load_model.hpp
  model_classes/classical.hpp
  model_classes/DBIntClass.hpp
  model_classes/exdc1.hpp
  model_classes/GainBlockClass.hpp
  model_classes/gensal.hpp
  model_classes/wsieg1.hpp
  model_classes/genrou.hpp
  model_classes/esst4b.hpp
  model_classes/esst1a.hpp
  model_classes/wshygp.hpp
  model_classes/ggov1.hpp
  model_classes/lvshbl.hpp
  model_classes/frqtpat.hpp
  model_classes/distr1.hpp
  model_classes/acmotor.hpp
  model_classes/ieel.hpp
  model_classes/motorw.hpp
  model_classes/psssim.hpp
  DESTINATION include/gridpack/applications/modules/dynamic_simulation_full_y
)

install(FILES 
  base_classes/base_generator_model.hpp
  base_classes/base_exciter_model.hpp
  base_classes/base_pss_model.hpp
  base_classes/base_governor_model.hpp
  base_classes/base_relay_model.hpp
  base_classes/base_load_model.hpp
  DESTINATION include/gridpack/applications/modules/dynamic_simulation_full_y/base_classes
)

install(FILES 
  model_classes/classical.hpp
  model_classes/DBIntClass.hpp
  model_classes/exdc1.hpp
  model_classes/GainBlockClass.hpp
  model_classes/gensal.hpp
  model_classes/wsieg1.hpp
  model_classes/genrou.hpp
  model_classes/esst4b.hpp
  model_classes/esst1a.hpp
  model_classes/wshygp.hpp
  model_classes/ggov1.hpp
  model_classes/lvshbl.hpp
  model_classes/frqtpat.hpp
  model_classes/distr1.hpp
  model_classes/acmotor.hpp
  model_classes/ieel.hpp
  model_classes/motorw.hpp
  model_classes/psssim.hpp
  DESTINATION include/gridpack/applications/modules/dynamic_simulation_full_y/model_classes
)

install(TARGETS 
  gridpack_dynamic_simulation_full_y_module
  DESTINATION lib
)
