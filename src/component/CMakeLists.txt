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
# Last Change: 2019-08-16 13:46:26 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_parallel
    gridpack_environment
    gridpack_math
    ${GA_LIBRARIES}
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
    include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_library(gridpack_components
  base_component.cpp
  data_collection.cpp
  optimization_ifc.cpp
  )
gridpack_set_library_version(gridpack_components)
add_dependencies(gridpack_components external_build)



# -------------------------------------------------------------
# component serialization tests
# -------------------------------------------------------------
add_executable(component_serialization_test test/component_serialization.cpp)
target_link_libraries(component_serialization_test gridpack_components ${target_libraries})
gridpack_add_unit_test(component_serialization component_serialization_test)

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  base_component.hpp  
  data_collection.hpp
  optimization_ifc.hpp
  DESTINATION include/gridpack/component
)

install(TARGETS 
  gridpack_components
  DESTINATION lib
)


