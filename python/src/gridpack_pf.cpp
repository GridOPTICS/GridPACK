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
    gpm.def_submodule("powerflow", "GridPACK power flow module");

  py::class_<gpf::PFAppModule> pfapp(pfm, "Powerflow");
  
  pfapp.doc() = ("Power flow application module");

  // minimal interface
  pfapp
    .def(py::init<>())
    .def("readNetwork",
         [](gpf::PFAppModule& self, gpu::Configuration& config, const int& idx) {

           gpp::Communicator world;
           boost::shared_ptr<gpf::PFNetwork> pfnet( new gpf::PFNetwork(world) );

           self.readNetwork(pfnet, &config, idx);
         }, "Read network specified in the configuration")
    .def("initialize", &gpf::PFAppModule::initialize,
         "Initialize power network")
    .def("reload", &gpf::PFAppModule::reload,
         "Reload network state and parameters from data collections")
    .def("nl_solve", &gpf::PFAppModule::nl_solve,
         "Solve the power flow problem using the math library non-linear solver")
    .def("solve", &gpf::PFAppModule::solve,
         "Solve the power flow problem using a custom Newton-Raphson loop")
    .def("saveData", &gpf::PFAppModule::saveData,
         "Save results of powerflow calculation to network data collection objects")
    .def("saveDataAlsotoOrg", &gpf::PFAppModule::saveDataAlsotoOrg,
         "Save results of powerflow calculation to data collection objects "
	 "also modify the original bus mag, ang, and the original generator "
         "PG QG in the datacollection")
    .def("exportPSSE23", &gpf::PFAppModule::exportPSSE23,
         "Export final configuration to PSS/E v23 formatted file")
    .def("exportPSSE33", &gpf::PFAppModule::exportPSSE33,
         "Export final configuration to PSS/E v33 formatted file")
    .def("exportPSSE34", &gpf::PFAppModule::exportPSSE34,
         "Export final configuration to PSS/E v34 formatted file")
    .def("suppressOutput", &gpf::PFAppModule::suppressOutput,
         "Suppress all output from power flow module")
    ;

  // Write results to standard output, or file
  pfapp
    .def("open",
         [](gpf::PFAppModule& self, const std::string& filename) {
           self.open(filename.c_str());
         },
         "Redirect power flow module output from standard out to a file")
    .def("close", &gpf::PFAppModule::close,
         "Close redirect output file")
    .def("write", &gpf::PFAppModule::write,
         "Write out results of powerflow calculation to standard output")
    ;
}
