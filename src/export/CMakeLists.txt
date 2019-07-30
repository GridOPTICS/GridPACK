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
# Created June 7, 2019 by Bruce Palmer
# Last Change: 2015-06-25 13:10:19 d3g096
# -------------------------------------------------------------

set(target_libraries
    ${GA_LIBRARIES}
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

if (PARMETIS_FOUND)
  include_directories(BEFORE ${PARMETIS_INCLUDE_DIR})
  set(target_libraries 
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${target_libraries}
    )
endif()

set(target_libraries
    gridpack_components
    gridpack_parallel
    gridpack_timer
    ${target_libraries})

include_directories(AFTER data_blocks)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  base_export.hpp
  PSSE33Export.hpp
  DESTINATION include/gridpack/export
)
install(FILES 
  data_blocks/export_2term33.hpp
  data_blocks/export_area33.hpp
  data_blocks/export_bus33.hpp
  data_blocks/export_facts33.hpp
  data_blocks/export_fxshnt33.hpp
  data_blocks/export_gen33.hpp
  data_blocks/export_iarea33.hpp
  data_blocks/export_icorr33.hpp
  data_blocks/export_line33.hpp
  data_blocks/export_load33.hpp
  data_blocks/export_msect33.hpp
  data_blocks/export_mterm33.hpp
  data_blocks/export_owner33.hpp
  data_blocks/export_swshnt33.hpp
  data_blocks/export_vscline33.hpp
  data_blocks/export_xform33.hpp
  data_blocks/export_zone33.hpp
  DESTINATION include/gridpack/export/data_blocks
)
