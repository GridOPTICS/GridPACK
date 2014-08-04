# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May 12, 2014 by William A. Perkins
# Last Change: 2014-05-15 09:14:52 d3g096
# -------------------------------------------------------------

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(BEFORE ${GRIDPACK_INCLUDE_DIRS})

add_library(gridpackf_network OBJECT
  network_c.cpp
  network_f.F90
)

#add_executable(hello_parallel
#  hello.F90
#  $<TARGET_OBJECTS:gridpackf_component>
#)

#target_link_libraries(hello_parallel
#  ${GRIDPACK_LIBS}
#)

