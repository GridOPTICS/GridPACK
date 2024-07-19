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
# Last Change: 2022-10-05 08:18:54 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_components
    gridpack_partition
    gridpack_parallel
    gridpack_stream
    gridpack_environment
    gridpack_math
    gridpack_configuration
    gridpack_timer
    gridpack_vsmatrix_components
    gridpack_ymatrix_components
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${PETSC_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_library(gridpack_voltage_stability_module
  vs_app_module.cpp
  vs_factory_module.cpp
  )

gridpack_set_library_version(gridpack_voltage_stability_module)

target_link_libraries(gridpack_voltage_stability_module
  gridpack_vsmatrix_components
  gridpack_partition
  gridpack_stream
  )


# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  vs_app_module.hpp
  vs_factory_module.hpp
  DESTINATION include/gridpack/applications/modules/voltage_stability
)

install(TARGETS 
  gridpack_voltage_stability_module
  DESTINATION lib
)
