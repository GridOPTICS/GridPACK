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
  gridpack_components
  gridpack_factory
  gridpack_partition
  gridpack_math
  gridpack_timer
  gridpack_expression
  gridpack_parallel
  gridpack_configuration
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_definitions(
  -DIL_STD
)

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

set(expression_src 
  variable.cpp
  expression.cpp
  optimizer.cpp
  lpfile_optimizer_implementation.cpp
)

if (USE_CPLEX)
  list(APPEND expression_src cplex_optimizer_implementation.cpp)
  add_definitions(-DHAVE_CPLEX)
  include_directories(AFTER ${CPLEX_INCLUDE_DIRS})
  list(APPEND target_libraries ${CPLEX_LIBRARIES})
elseif(USE_GLPK)
  list(APPEND expression_src glpk_optimizer_implementation.cpp)
  add_definitions(-DHAVE_GLPK)
  include_directories(BEFORE ${GLPK_INCLUDE_DIR})
  list(APPEND target_libraries ${GLPK_LIBRARY})
elseif(USE_GUROBI)
  include_directories(AFTER ${GUROBI_INCLUDE_DIRS})
  list(APPEND target_libraries ${GUROBI_LIBRARIES})
else()
  message(WARNING "Optimization library not specified")
endif()

add_library(gridpack_expression
  ${expression_src}
)

list(APPEND target_libraries
  ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
  ${GA_LIBRARIES} 
  ${Boost_LIBRARIES} 
  ${PETSC_LIBRARIES}
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
# optimizer test suite
# -------------------------------------------------------------
add_executable(optimizer_test test/optimizer_test.cpp)
target_link_libraries(optimizer_test ${target_libraries})

gridpack_add_unit_test("optimizer" optimizer_test)

# -------------------------------------------------------------
# TEST: 
# A simple optimization problem (3 generators)
# -------------------------------------------------------------
if(USE_CPLEX)
 add_executable(test_optim 
   test/test_optim.cpp
 )
 target_link_libraries(test_optim ${target_libraries})

 add_executable(test_expression 
   test/test_expression.cpp
 )
 target_link_libraries(test_expression ${target_libraries})
endif()


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
   
   COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_SOURCE_DIR}/test_app/loads.csv
   ${CMAKE_CURRENT_BINARY_DIR} 

   DEPENDS
   ${CMAKE_CURRENT_SOURCE_DIR}/test_app/loads.csv
 )
 add_dependencies(uc.x uc.x.input)

if(USE_GUROBI)
 add_executable(uc_gurobi.x
   test_gurobi/uc_app.cpp
   test_gurobi/uc_components.cpp
   test_gurobi/uc_main.cpp
 )

 target_link_libraries(uc_gurobi.x ${target_libraries})

 add_custom_target(uc_gurobi.x.input

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
 add_dependencies(uc_gurobi.x uc_gurobi.x.input)
endif()


# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  variable.hpp  
  expression.hpp
  optimization.hpp
  DESTINATION include/gridpack/optimization
)

install(TARGETS 
  gridpack_expression
  DESTINATION lib
)
