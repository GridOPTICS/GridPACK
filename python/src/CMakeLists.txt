# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
#
# Copyright (c) 2013 Battelle Memorial Institute
# Licensed under modified BSD License. A copy of this license can be found
# in the LICENSE file in the top level directory of this distribution.
#
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created February 17, 2020 by Perkins
# Last Change: 2025-08-11 12:30:14 d3g096
# Added gridpack_cli.py + MPI include    
# Updated in March 2026, by Yousu Chen
# -------------------------------------------------------------

pybind11_add_module(_gridpack
  gridpack.cpp
  gridpack_ds.cpp
  gridpack_emt.cpp
  gridpack_pf.cpp
  gridpack_se.cpp
  )

set(parallel_scripts
  ${CMAKE_CURRENT_SOURCE_DIR}/hello.py
  ${CMAKE_CURRENT_SOURCE_DIR}/task_manager.py
  ${CMAKE_CURRENT_SOURCE_DIR}/hadrec.py
  ${CMAKE_CURRENT_SOURCE_DIR}/emt.py
  ${CMAKE_CURRENT_SOURCE_DIR}/pf.py
  ${CMAKE_CURRENT_SOURCE_DIR}/stes.py
)

if (ENABLE_ENVIRONMENT_FROM_COMM)
  list(APPEND parallel_scripts ${CMAKE_CURRENT_SOURCE_DIR}/hello_comm.py)
endif()

add_definitions(${GRIDPACK_DEFINITIONS})
target_include_directories(_gridpack PRIVATE ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(_gridpack PRIVATE ${GRIDPACK_LIBS})

install(FILES
  ${parallel_scripts}
  DESTINATION bin
  )


