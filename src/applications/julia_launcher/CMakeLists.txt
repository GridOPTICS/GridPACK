#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be found
#     in the LICENSE file in the top level directory of this distribution.
#
# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# Julia Launcher Application
# Distributes OPF tasks to Julia workers using GridPACK's TaskManager
# -------------------------------------------------------------
# Created 2026-01-09
# -------------------------------------------------------------

# This application only needs parallel components (no power flow modules)
set(target_libraries
    gridpack_environment
    gridpack_parallel
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# Build the julia_launcher executable
add_executable(julia_launcher.x
   jl_driver.cpp
   jl_main.cpp
)

target_link_libraries(julia_launcher.x ${target_libraries})

# Install the executable
install(TARGETS julia_launcher.x DESTINATION bin)

# Install source files as example
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/jl_driver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/jl_driver.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/jl_main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
  DESTINATION share/gridpack/example/julia_launcher
)
