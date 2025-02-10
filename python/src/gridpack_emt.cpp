// -------------------------------------------------------------
// file: gridpack_emt.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created February 10, 2025 by Perkins
// -------------------------------------------------------------

#include "common.hpp"
#include "gridpack/applications/modules/emt/emt.hpp"

// -------------------------------------------------------------
// init_gridpack_emt
// Inter facet to Emt application module
// -------------------------------------------------------------
void
init_gridpack_emt(py::module& gpm)
{
  py::module emtm =
    gpm.def_submodule("emt", "GridPACK EMT Application module");

  py::class_<Emt> emtapp(emtm, "EMT");
  emtapp
    .def(py::init<>())
    .def(py::init<gpp::Communicator&>())
    .def("rank", &Emt::rank)
    .def("size", &Emt::size)
    .def("setup", &Emt::setup)
    .def("initialize", &Emt::initialize)
    .def("solve", &Emt::solve)
    .def("readnetworkdatafromconfig", &Emt::readnetworkdatafromconfig)
    .def("solvepowerflow", &Emt::solvepowerflow)
    .def("setconfigurationfile",
         [](Emt& self, const std::string& s) {
           self.setconfigurationfile(s.c_str());
         })
    ;

}  

