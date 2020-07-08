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
# Last Change: 2015-06-25 13:10:19 d3g096
# -------------------------------------------------------------

set(target_libraries
    ${GA_LIBRARIES}
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

if (PARMETIS_FOUND)
  include_directories(BEFORE ${PARMETIS_INCLUDE_DIR})
  set(target_libraries 
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${target_libraries}
    )
endif()

set(target_libraries
    gridpack_environment
    gridpack_math
    gridpack_partition
    gridpack_components
    gridpack_parallel
    gridpack_timer
    gridpack_stream
    ${target_libraries})

include_directories(AFTER parser_classes)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()
# -------------------------------------------------------------
# TEST: parser_test
# A simple program to check input functionality
# -------------------------------------------------------------
#add_executable(PTI23_test test/PTI23_test.cpp)
#target_link_libraries(PTI23_test ${target_libraries})

add_executable(parser_test test/parser_test.cpp)
target_link_libraries(parser_test ${target_libraries})

add_custom_target(test_parser_input

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/test/test.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/test/parser_data.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GridPACK_SOURCE_DIR}/applications/data_sets/raw/IEEE14.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/test/table.dat
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/test/test.xml
  ${CMAKE_CURRENT_SOURCE_DIR}/test/parser_data.raw
  ${GridPACK_SOURCE_DIR}/applications/data_sets/raw/IEEE14.raw
  ${CMAKE_CURRENT_SOURCE_DIR}/test/table.dat
)
add_dependencies(parser_test test_parser_input)
#add_dependencies(PTI23_test test_parser_input)

gridpack_add_run_test(parser_test parser_test "IEEE14.raw")

# -------------------------------------------------------------
# TEST: hash_distr_test
# -------------------------------------------------------------
add_executable(hash_distr_test test/hash_distr_test.cpp)
target_link_libraries(hash_distr_test ${target_libraries})

gridpack_add_unit_test(hash_distr_test hash_distr_test)

# -------------------------------------------------------------
# TEST: bus_table_test
# -------------------------------------------------------------
add_executable(bus_table_test test/bus_table_test.cpp)
target_link_libraries(bus_table_test ${target_libraries})

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  dictionary.hpp
  PTI23_parser.hpp
  PTI33_parser.hpp
  GOSS_parser.hpp
  MAT_parser.hpp
  hash_distr.hpp
  base_parser.hpp
  base_pti_parser.hpp
  bus_table.hpp
  DESTINATION include/gridpack/parser
)
install(FILES 
  parser_classes/gencls.hpp
  parser_classes/gensal.hpp
  parser_classes/genrou.hpp
  parser_classes/wsieg1.hpp
  parser_classes/exdc1.hpp
  parser_classes/esst1a.hpp
  parser_classes/esst4b.hpp
  parser_classes/ggov1.hpp
  parser_classes/tgov1.hpp
  parser_classes/wshygp.hpp
  parser_classes/lvshbl.hpp
  parser_classes/frqtpat.hpp
  parser_classes/distr1.hpp
  parser_classes/cim6bl.hpp
  parser_classes/acmtblu1.hpp
  parser_classes/ieelbl.hpp
  parser_classes/cmldblu1.hpp
  parser_classes/psssim.hpp
  DESTINATION include/gridpack/parser/parser_classes
)
