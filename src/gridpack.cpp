// -------------------------------------------------------------
// file: gridpack.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created January 24, 2020 by Perkins
// Last Change: 2020-01-24 13:01:58 d3g096
// -------------------------------------------------------------

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <gridpack/environment/environment.hpp>
#include <gridpack/parallel/communicator.hpp>
namespace gp = gridpack;
namespace gpp = gridpack::parallel;

PYBIND11_MODULE(gridpack, gpm) {
  gpm.doc() = "GridPACK module";
}
