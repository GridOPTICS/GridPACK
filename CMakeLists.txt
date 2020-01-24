cmake_minimum_required(VERSION 2.8.12)
project(gridpack-hadrec)

find_package(pybind11 REQUIRED)

if (NOT GRIDPACK_DIR)
  set(GRIDPACK_DIR "/Users/d3g096/Projects/GridPACK/gridpack-hadrec"
    CACHE PATH "GridPACK installation directory")
endif()

include("${GRIDPACK_DIR}/lib/GridPACK.cmake")

gridpack_setup()

add_definitions(${GRIDPACK_DEFINITIONS})
include_directories(BEFORE ${GRIDPACK_INCLUDE_DIRS})


add_subdirectory(src)
