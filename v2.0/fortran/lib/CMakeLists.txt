# -*- mode: cmake -*-------------------------------------------
#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can
#     be found in the LICENSE file in the top level directory of this
#     distribution.
#
# -------------------------------------------------------------
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May 19, 2014 by William A. Perkins
# Last Change: 2014-05-19 09:50:02 d3g096
# -------------------------------------------------------------

add_library(gridpackf
  $<TARGET_OBJECTS:gridpackf_parallel>
  $<TARGET_OBJECTS:gridpackf_math>
  $<TARGET_OBJECTS:gridpackf_network>
  $<TARGET_OBJECTS:gridpackf_component>
)

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------

install(TARGETS 
  gridpackf
  DESTINATION lib
)
