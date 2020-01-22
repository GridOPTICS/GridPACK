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
# Last Change: 2019-12-16 12:57:42 d3g096
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
    gridpack_pfmatrix_components
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

add_library(gridpack_powerflow_module
  pf_app_module.cpp
  pf_factory_module.cpp
  )

gridpack_set_library_version(gridpack_powerflow_module)

target_link_libraries(gridpack_powerflow_module
  gridpack_pfmatrix_components
  gridpack_partition
  gridpack_math
  gridpack_stream
  )


# -------------------------------------------------------------
# component serialization tests
# -------------------------------------------------------------
# target_link_libraries(gridpack_powerflow_module ${target_libraries})

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  pf_app_module.hpp
  pf_factory_module.hpp
  DESTINATION include/gridpack/applications/modules/powerflow
)

install(TARGETS 
  gridpack_powerflow_module
  DESTINATION lib
)
