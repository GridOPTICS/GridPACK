# -*- mode: cmake -*-
#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be
#     found
#     in the LICENSE file in the top level directory of this distribution.
#
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created June 11, 2013 by William A. Perkins
# Last Change: 2015-09-01 14:25:11 d3g096
# -------------------------------------------------------------

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(gridpack_expression
  variable.cpp
  expression.cpp
  optimizer.cpp
  cplex_optimizer_implementation.cpp
)

set(target_libraries 
  gridpack_expression
  gridpack_configuration 
  gridpack_parallel 
  ${GA_LIBRARIES}
  ${Boost_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
)

add_executable(variable_test variable_test.cpp)
target_link_libraries(variable_test ${target_libraries})

add_executable(expression_test expression_test.cpp)
target_link_libraries(expression_test ${target_libraries})

add_executable(optimizer_test optimizer_test.cpp)
target_link_libraries(optimizer_test ${target_libraries})

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  variable.hpp  
  expression.hpp
  DESTINATION include/gridpack/expression
)

install(TARGETS 
  gridpack_expression
  DESTINATION lib
)
