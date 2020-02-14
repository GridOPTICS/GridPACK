// -------------------------------------------------------------
// file: test_mpi_wrap.cpp
// -------------------------------------------------------------
/*
 * Copyright (c) 2013 Battelle Memorial Institute
 * Licensed under modified BSD License. A copy of this license can be found
 * in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// Created February 14, 2020 by Perkins
// Last Change: 2020-02-14 07:23:45 d3g096
// -------------------------------------------------------------


#include <vector>

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "test_mpi.hpp"

PYBIND11_MODULE(Test, m) {
  py::class_<TestSerialMPI>(m, "SerialMPI")
    .def(py::init<>())
    ;
}

