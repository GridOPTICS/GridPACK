# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May  6, 2013 by William A. Perkins
# Last Change: 2013-08-09 11:41:42 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_components
    gridpack_factory
    gridpack_partition
    gridpack_parallel
    gridpack_math
    gridpack_configuration
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${PETSC_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_library(powerflow_app
  pf_components.cpp
  pf_factory.cpp
)

add_executable(powerflow.x
   pf_app.cpp
   pf_components.cpp
   pf_factory.cpp
   pf_main.cpp
)

target_link_libraries(powerflow.x ${target_libraries})


# -------------------------------------------------------------
# TEST: pf_network_test
# -------------------------------------------------------------
add_executable(pf_network_test test/pf_network_test.cpp)
target_link_libraries(pf_network_test powerflow_app ${target_libraries})
gridpack_add_unit_test(pf_network_test pf_network_test)

