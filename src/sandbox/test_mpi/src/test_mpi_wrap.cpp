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
// Last Change: 2020-02-14 14:34:22 d3g096
// -------------------------------------------------------------


#include <pybind11/pybind11.h>
namespace py = pybind11;
#include <pybind11/stl.h>

#include "test_mpi.hpp"

// Some temporary hacks

// #define RHEL_OPENMPI_HACK 1
#ifdef RHEL_OPENMPI_HACK

// This stupidity is needed on RHEL7 with stock OpenMPI packages
// installed.

#include <dlfcn.h>

// -------------------------------------------------------------
// stupid_openmpi_hack
// from https://github.com/baidu-research/tensorflow-allreduce/issues/4
// -------------------------------------------------------------
static void 
stupid_openmpi_hack(void)
{
  void *handle = NULL;
  int mode = RTLD_NOW | RTLD_GLOBAL;

  // GNU/Linux and others 
#ifdef RTLD_NOLOAD
      mode |= RTLD_NOLOAD;
#endif
  if (!handle) handle = dlopen("libmpi.so.20", mode);
  if (!handle) handle = dlopen("libmpi.so.12", mode);
  if (!handle) handle = dlopen("libmpi.so.1", mode);
  if (!handle) handle = dlopen("libmpi.so.0", mode);
  if (!handle) handle = dlopen("libmpi.so", mode);
}

#endif

PYBIND11_MODULE(Test, m) {
#ifdef RHEL_OPENMPI_HACK
  stupid_openmpi_hack();
#endif
  py::class_<TestSerialMPI>(m, "SerialMPI")
    .def(py::init<>())
    .def("show", &TestSerialMPI::show)
    ;
}

