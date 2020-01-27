
pybind11_add_module(gridpack gridpack.cpp)

add_custom_target(parallel_scripts
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/hello.py ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS 
  ${CMAKE_CURRENT_SOURCE_DIR}/hello.py
  )

add_dependencies(gridpack parallel_scripts)

