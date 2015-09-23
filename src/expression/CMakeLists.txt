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
# Last Change: 2015-09-16 10:19:43 d3g096
# -------------------------------------------------------------

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(expression_src 
  variable.cpp
  expression.cpp
  optimizer.cpp
  lpfile_optimizer_implementation.cpp
)

set(target_libraries 
  gridpack_expression
  gridpack_configuration 
  gridpack_parallel 
)

if(USE_CPLEX)
  list(APPEND expression_src cplex_optimizer_implementation.cpp)
  add_definitions(-DHAVE_CPLEX)
elseif(USE_GLPK)
  list(APPEND expression_src glpk_optimizer_implementation.cpp)
  add_definitions(-DHAVE_GLPK)
  include_directories(BEFORE ${GLPK_INCLUDE_DIR})
  list(APPEND target_libraries ${GLPK_LIBRARY})
else()
  message(WARNING "Optimization library not specified")
endif()

add_library(gridpack_expression
  ${expression_src}
)

list(APPEND target_libraries
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
