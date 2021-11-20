cmake_minimum_required(VERSION 2.8.12)
project(gridpack-hadrec)

if (DEFINED ENV{RHEL_OPENMPI_HACK})
  add_definitions (-DRHEL_OPENMPI_HACK=1)
endif()
  
# -------------------------------------------------------------
# PyBind11
# We can use a (correctly) installed pybind11 or use source tree in the
# top directory.  Get the pybind11 source using:
#
# git clone -b v2.4 https://github.com/pybind/pybind11.git
# -------------------------------------------------------------
if (EXISTS ${PROJECT_SOURCE_DIR}/pybind11)
  add_subdirectory(pybind11)
else()
  find_package(pybind11 REQUIRED)
endif()

# -------------------------------------------------------------
# GridPACK
# -------------------------------------------------------------
if (DEFINED ENV{GRIDPACK_DIR})
  set(GRIDPACK_DIR "$ENV{GRIDPACK_DIR}"
    CACHE PATH "GridPACK installation directory")
endif()

if (NOT GRIDPACK_DIR)
  set(GRIDPACK_DIR "/Users/d3g096/Projects/GridPACK/gridpack-hadrec"
    CACHE PATH "GridPACK installation directory")
endif()

include("${GRIDPACK_DIR}/lib/GridPACK.cmake")

gridpack_setup()

add_definitions(${GRIDPACK_DEFINITIONS})
include_directories(BEFORE ${GRIDPACK_INCLUDE_DIRS})


add_subdirectory(src)
