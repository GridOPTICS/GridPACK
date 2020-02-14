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
# -------------------------------------------------------------
find_package(pybind11 REQUIRED)

add_subdirectory(src)
