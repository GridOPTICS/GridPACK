# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created February 17, 2020 by Perkins
# Last Change: 2022-11-03 10:41:21 d3g096
# -------------------------------------------------------------

pybind11_add_module(gridpack gridpack.cpp)

add_custom_target(parallel_scripts
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/hello.py ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/hello_comm.py ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/task_manager.py ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/hadrec.py ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS 
  ${CMAKE_CURRENT_SOURCE_DIR}/hello.py
  ${CMAKE_CURRENT_SOURCE_DIR}/hello_comm.py
  ${CMAKE_CURRENT_SOURCE_DIR}/task_manager.py
  ${CMAKE_CURRENT_SOURCE_DIR}/hadrec.py
  )

add_definitions(${GRIDPACK_DEFINITIONS})
add_dependencies(gridpack parallel_scripts)
target_link_libraries(gridpack PRIVATE ${GRIDPACK_LIBS})

