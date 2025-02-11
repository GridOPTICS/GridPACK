// -------------------------------------------------------------
// file: gridpack_pf.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
/*
 * Copyright (c) 2013 Battelle Memorial Institute
 * Licensed under modified BSD License. A copy of this license can be found
 * in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created February 10, 2025 by Perkins
// -------------------------------------------------------------

#include "common.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"

namespace gpf = gridpack::powerflow;

// -------------------------------------------------------------
// init_gridpack_pf
// Interface to powerflow application module
// -------------------------------------------------------------
void
init_gridpack_pf(py::module& gpm)
{
  py::module pfm =
    gpm.def_submodule("powerflow", "GridPACK Powerflow Application module");

  py::class_<gpf::PFAppModule> pfapp(pfm, "powerflow");


  // minimal interface
  pfapp
    .def(py::init<>())
    .def("readNetwork",
         [](gpf::PFAppModule& self, gpu::Configuration& config, const int& idx) {

           gpp::Communicator world;
           boost::shared_ptr<gpf::PFNetwork> pfnet( new gpf::PFNetwork(world) );

           self.readNetwork(pfnet, &config, idx);
         })
    .def("initialize", &gpf::PFAppModule::initialize)
    .def("nl_solve", &gpf::PFAppModule::nl_solve)
    .def("solve", &gpf::PFAppModule::solve)
    .def("write", &gpf::PFAppModule::write)
    .def("saveData", &gpf::PFAppModule::saveData)
    .def("exportPSSE23", &gpf::PFAppModule::exportPSSE23)
    .def("exportPSSE33", &gpf::PFAppModule::exportPSSE33)
    .def("exportPSSE34", &gpf::PFAppModule::exportPSSE34)
    ;
}
