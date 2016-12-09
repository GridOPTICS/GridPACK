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

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_definitions(
  -DIL_STD
)

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

set(expression_src 
  optimizer.cpp
  file_optimizer_implementation.cpp
  lpfile_optimizer_implementation.cpp
  julia_optimizer_implementation.cpp
)

set(optimizer_libraries)

if (USE_CPLEX)
  list(APPEND expression_src cplex_optimizer_implementation.cpp)
  add_definitions(-DHAVE_CPLEX)
  include_directories(AFTER ${CPLEX_INCLUDE_DIRS})
  list(APPEND optimizer_libraries ${CPLEX_LIBRARIES})
endif()
if (USE_GLPK)
  list(APPEND expression_src glpk_optimizer_implementation.cpp)
  add_definitions(-DHAVE_GLPK)
  include_directories(BEFORE ${GLPK_INCLUDE_DIR})
  list(APPEND optimizer_libraries ${GLPK_LIBRARY})
endif()
if (USE_GUROBI)
  include_directories(AFTER ${GUROBI_INCLUDE_DIRS})
  list(APPEND optimizer_libraries ${GUROBI_LIBRARIES})
endif()

add_library(gridpack_optimum
  ${expression_src}
)


# target_link_libraries(gridpack_optimum
#   gridpack_parallel
#   gridpack_configuration
#   ${optimizer_libraries}
#   ${GA_LIBRARIES} 
#   ${Boost_LIBRARIES} 
#   ${PETSC_LIBRARIES}
#   ${MPI_CXX_LIBRARIES}
# )
  
set(target_libraries
  gridpack_optimum
  gridpack_components
  gridpack_partition
  gridpack_math
  gridpack_timer
  gridpack_expression
  gridpack_parallel
  gridpack_configuration
  ${optimizer_libraries}
  ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
  ${GA_LIBRARIES} 
  ${Boost_LIBRARIES} 
  ${PETSC_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
)

# -------------------------------------------------------------
# optimizer test suite
# -------------------------------------------------------------
add_executable(optimizer_test test/optimizer_test.cpp)
target_link_libraries(optimizer_test ${target_libraries})

gridpack_add_unit_test("optimizer" optimizer_test)

add_custom_target(optimizer_test_input
  
  COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/test/gridpack.xml
  ${CMAKE_CURRENT_BINARY_DIR} 
  
  DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/test/gridpack.xml
  )
add_dependencies(optimizer_test optimizer_test_input)

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
  optimization.hpp
  DESTINATION include/gridpack/optimization
)

install(TARGETS 
  gridpack_optimum
  DESTINATION lib
)
