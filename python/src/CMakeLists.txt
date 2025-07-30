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
# Last Change: 2025-03-11 14:04:06 d3g096
# -------------------------------------------------------------

pybind11_add_module(gridpack
  gridpack.cpp
  gridpack_ds.cpp
  gridpack_emt.cpp
  gridpack_pf.cpp
  )

set(parallel_scripts 
  ${CMAKE_CURRENT_SOURCE_DIR}/hello.py
  ${CMAKE_CURRENT_SOURCE_DIR}/task_manager.py
  ${CMAKE_CURRENT_SOURCE_DIR}/hadrec.py
  ${CMAKE_CURRENT_SOURCE_DIR}/emt.py
  ${CMAKE_CURRENT_SOURCE_DIR}/pf.py
)

if (ENABLE_ENVIRONMENT_FROM_COMM)
  list(APPEND parallel_scripts ${CMAKE_CURRENT_SOURCE_DIR}/hello_comm.py)
endif()

add_definitions(${GRIDPACK_DEFINITIONS})
target_link_libraries(gridpack PRIVATE ${GRIDPACK_LIBS})

install(FILES
  ${parallel_scripts}
  DESTINATION bin
  )


