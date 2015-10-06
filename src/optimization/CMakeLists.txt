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
# Last Change: 2015-07-30 13:20:39 d3g096
# -------------------------------------------------------------

set(target_libraries
  gridpack_components
  gridpack_factory
  gridpack_partition
  gridpack_parallel
  gridpack_math
  gridpack_timer
  ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
  ${GA_LIBRARIES} 
  ${Boost_LIBRARIES} 
  ${MPI_CXX_LIBRARIES}
  ${CPLEX_LIBRARIES}
  ${PETSC_LIBRARIES}
)

add_definitions(
     -DIL_STD
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif(
)
if (CPLEX_FOUND)
  include_directories(AFTER ${CPLEX_INCLUDE_DIRS})
endif()
# -------------------------------------------------------------
# TEST: greetings
# A simple MPI program
# -------------------------------------------------------------
add_executable(test_expression
   test/test_expression.cpp
)
target_link_libraries(test_expression ${target_libraries})

add_executable(uc.x
   test_app/uc_app.cpp
   test_app/uc_components.cpp
   test_app/uc_main.cpp
)

target_link_libraries(uc.x ${target_libraries})

add_custom_target(uc.x.input

   COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_SOURCE_DIR}/test_app/uc_test.raw
   ${CMAKE_CURRENT_BINARY_DIR} 

   DEPENDS
   ${CMAKE_CURRENT_SOURCE_DIR}/test_app/uc_test.raw

   COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_SOURCE_DIR}/test_app/gen.uc
   ${CMAKE_CURRENT_BINARY_DIR} 

   DEPENDS
   ${CMAKE_CURRENT_SOURCE_DIR}/test_app/gen.uc

   COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_SOURCE_DIR}/test_app/loads.txt
   ${CMAKE_CURRENT_BINARY_DIR} 

   DEPENDS
   ${CMAKE_CURRENT_SOURCE_DIR}/test_app/loads.txt
)
add_dependencies(uc.x uc.x.input)

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  optimization.hpp
  DESTINATION include/gridpack/optimization
)

