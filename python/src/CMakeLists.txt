# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created February 17, 2020 by Perkins
# Last Change: 2023-12-12 10:26:16 d3g096
# -------------------------------------------------------------

pybind11_add_module(gridpack gridpack.cpp)

set(parallel_scripts 
  ${CMAKE_CURRENT_SOURCE_DIR}/hello.py
  ${CMAKE_CURRENT_SOURCE_DIR}/task_manager.py
  ${CMAKE_CURRENT_SOURCE_DIR}/hadrec.py
)

if (ENABLE_ENVIRONMENT_FROM_COMM)
  list(APPEND parallel_scripts ${CMAKE_CURRENT_SOURCE_DIR}/hello_comm.py)
endif()

add_definitions(${GRIDPACK_DEFINITIONS})
target_link_libraries(gridpack PRIVATE ${GRIDPACK_LIBS})

install(FILES
  ${parallel_scripts}
  DESTINATION bin)


