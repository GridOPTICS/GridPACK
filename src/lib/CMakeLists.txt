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
# Created October 17, 2013 by Bruce Palmer
# Last Change: 2014-01-15 13:29:41 d3g096
# -------------------------------------------------------------

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(GRIDPACK_TP_INCLUDE_FLAGS)
  
set(GRIDPACK_TP_LIBS)


# -------------------------------------------------------------
#  Add a single library for all GridPACK functions
# -------------------------------------------------------------
# source files
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/component/base_component.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/component/data_collection.cpp)

set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/configuration/configuration.cpp)

set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/implementation_visitor.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/linear_solver.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/linear_solver_implementation.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/matrix.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/matrix_implementation.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/newton_raphson_solver.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/newton_raphson_solver_implementation.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/nonlinear_solver.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/nonlinear_solver_implementation.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/nonlinear_solver_interface.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/vector.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/vector_implementation.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/petsc/linear_solver.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/petsc/math.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/petsc/matrix.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/petsc/nonlinear_solver.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/petsc/petsc_configuration.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/petsc/petsc_exception.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/petsc/petsc_linear_solver_implementation.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/petsc/petsc_matrix_implementation.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/petsc/petsc_matrix_operations.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/petsc/petsc_nonlinear_solver_implementation.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/petsc/petsc_vector_implementation.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/math/petsc/vector.cpp)

set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/parallel/distributed.cpp)

set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/partition/adjacency_list.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/partition/graph_partitioner.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/partition/graph_partitioner_implementation.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/partition/simple_adjacency.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/partition/parmetis/graph_partitioner.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/partition/parmetis/parmetis_graph_partitioner_impl.cpp)
set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/partition/parmetis/parmetis_graph_wrapper.cpp)

set (SOURCES ${SOURCES} ${GRIDPACK_SRC_DIR}/timer/coarse_timer.cpp)

# header directories
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/component)
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/configuration)
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/factory)
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/mapper)
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/math)
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/math/petsc)
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/network)
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/parallel)
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/parser)
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/partition)
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/serial_io)
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/timer)
set(HEADERS ${HEADERS} ${GRIDPACK_SRC_DIR}/utilities)
#include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${HEADERS})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()
if (PARMETIS_FOUND)
  include_directories(BEFORE ${PARMETIS_INCLUDE_DIR})
endif()
if (PETSC_FOUND)
  include_directories(AFTER ${PETSC_INCLUDES})
  add_definitions (${PETSC_DEFINITIONS})
endif()

add_library(gridpack ${SOURCES})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GridPACK.cmake.in 
  ${CMAKE_CURRENT_BINARY_DIR}/GridPACK.cmake @ONLY)

install(FILES 
  ${CMAKE_CURRENT_BINARY_DIR}/GridPACK.cmake
  DESTINATION lib
)
