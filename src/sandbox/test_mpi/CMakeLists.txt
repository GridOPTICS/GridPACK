cmake_minimum_required(VERSION 2.8.12)
project(test CXX)


# -------------------------------------------------------------
# MPI compiler
# -------------------------------------------------------------
message(STATUS "Checking MPI ...")
find_package (MPI REQUIRED)
include_directories(BEFORE ${MPI_CXX_INCLUDE_PATH})

message(STATUS "MPI_CXX_INCLUDE_PATH: ${MPI_CXX_INCLUDE_PATH}")
message(STATUS "MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")


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
add_subdirectory(src)
