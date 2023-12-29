cmake_minimum_required(VERSION 2.8.12)
project(gridpack-hadrec)

if (DEFINED ENV{RHEL_OPENMPI_HACK})
  add_definitions (-DRHEL_OPENMPI_HACK=1)
endif()

# -------------------------------------------------------------
# Python
# -------------------------------------------------------------
find_package(Python REQUIRED COMPONENTS Interpreter Development)
message("*****************************************")
message("Python_ROOT " ${Python_ROOT})
find_package(Python REQUIRED COMPONENTS Interpreter Development)

message("Python_FOUND " ${Python_FOUND})
message("Python_Interpreter_FOUND " ${Python_Interpreter_FOUND})
message("Python_Development_FOUND " ${Python_Development_FOUND})
   
message("Python_LIBRARIES " ${Python_LIBRARIES})
message("Python_LIBRARY_DIRS " ${Python_LIBRARY_DIRS})
message("Python_INCLUDE_DIRS " ${Python_INCLUDE_DIRS})
message("Python_LINK_OPTIONS " ${Python_LINK_OPTIONS})
message("Python_EXECUTABLE " ${Python_EXECUTABLE})
message("Python_INTERPRETER_ID " ${Python_INTERPRETER_ID})
  
message("Python_VERSION " ${Python_VERSION})
message("Python_VERSION_MAJOR " ${Python_VERSION_MAJOR})
message("Python_VERSION_MINOR " ${Python_VERSION_MINOR})
message("*****************************************")

  
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
# mpi4py
# -------------------------------------------------------------
execute_process(
  COMMAND ${PYTHON_EXECUTABLE} -c "import mpi4py;print(mpi4py.get_include())"
  OUTPUT_VARIABLE MPI4PY_INCLUDE
  RESULT_VARIABLE EXIT_CODE
)

if (NOT ${EXIT_CODE} EQUAL 0)
  message(FATAL_ERROR "The mpi4py Python module is required but cannot be found")
endif()

message(STATUS "MPI4PY_INCLUDE=${MPI4PY_INCLUDE}")
include_directories(BEFORE ${MPI4PY_INCLUDE})

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
