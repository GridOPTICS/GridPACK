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
# Last Change: 2015-10-06 10:15:22 d3g096
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
add_definitions(
  -DIL_STD
)


if(USE_CPLEX)
  list(APPEND expression_src cplex_optimizer_implementation.cpp)
  add_definitions(-DHAVE_CPLEX)
  include_directories(AFTER ${CPLEX_INCLUDE_DIRS})
  list(APPEND target_libraries ${CPLEX_LIBRARIES})
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


# -------------------------------------------------------------
# variable test suite
# -------------------------------------------------------------
add_executable(variable_test test/variable_test.cpp)
target_link_libraries(variable_test ${target_libraries})

gridpack_add_unit_test("optimization_variables" variable_test)


# -------------------------------------------------------------
# expression test suite
# -------------------------------------------------------------
add_executable(expression_test expression_test.cpp)
target_link_libraries(expression_test ${target_libraries})

# -------------------------------------------------------------
# optimizer test suite
# -------------------------------------------------------------
add_executable(optimizer_test test/optimizer_test.cpp)
target_link_libraries(optimizer_test ${target_libraries})

gridpack_add_unit_test("optimizer" optimizer_test)

add_executable(uc_test uc_test.cpp)
target_link_libraries(uc_test ${target_libraries})

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
