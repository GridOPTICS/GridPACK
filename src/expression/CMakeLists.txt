# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be
#     found
#     in the LICENSE file in the top level directory of this distribution.
#
# -------------------------------------------------------------
# Created 
# Last Change:
# -------------------------------------------------------------
set(target_libraries
  gridpack_parallel
  ${Boost_LIBRARIES} 
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(expression_src 
  variable.cpp
  expression.cpp
  functions.cpp
)

add_library(gridpack_expression
  ${expression_src}
)

list(APPEND target_libraries
  gridpack_expression
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
add_executable(expression_test test/expression_test.cpp)
target_link_libraries(expression_test ${target_libraries})

gridpack_add_unit_test("optimization_expressions" variable_test)

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  variable.hpp  
  expression.hpp
  functions.hpp
  DESTINATION include/gridpack/expression
)

install(TARGETS 
  gridpack_expression
  DESTINATION lib
)
