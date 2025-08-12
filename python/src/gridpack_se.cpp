// -------------------------------------------------------------
// file: gridpack_se.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
/*
 * Copyright (c) 2013 Battelle Memorial Institute
 * Licensed under modified BSD License. A copy of this license can be found
 * in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August  6, 2025 by Perkins
// Last Change: 2025-08-11 13:00:19 d3g096
// -------------------------------------------------------------

#include "common.hpp"
#include "gridpack/applications/modules/state_estimation/se_app_module.hpp"

namespace gse = gridpack::state_estimation;

// -------------------------------------------------------------
// init_gridpack_se
// Interface to the GridPACK state estimation application module
// -------------------------------------------------------------
void
init_gridpack_se(py::module& gpm)
{
  py::module sem =
    gpm.def_submodule("state_estimation", "GridPACK state estimation module");

  py::class_<gse::SEAppModule> seapp(sem, "SEApp");
  
  seapp.doc() = ("State estimation application module");

  seapp
    .def(py::init<>())
    .def("readNetwork",
         [](gse::SEAppModule& self, gpu::Configuration& config) {
           gpp::Communicator world;
           boost::shared_ptr<gse::SENetwork> senet(new gse::SENetwork(world));
           self.readNetwork(senet, &config);
         }, R"eof(
Read in and partition the network. The input file is read directly
from the state_estimation block in the specified configuration.

Parameters:
    config (gridpack.Configuration): power flow problem configuration usually read from an XML file
)eof")

    .def("initialize", &gse::SEAppModule::initialize,
         R"eof( 
Set up exchange buffers and other internal parameters and initialize
network components using data from data collection
)eof")
    
    .def("readMeasurements", &gse::SEAppModule::readMeasurements,
         R"eof(
Read branch and bus measurements. These will come from a separate
file.  The name of this file comes from the input configuration from
which the network was read.
)eof")

    .def("solve", &gse::SEAppModule::solve,
         R"eof( 
Solve the state estimation problem
)eof")

    .def("write", &gse::SEAppModule::write,
         R"eof(
Write final results of state estimation calculation to standard output
)eof")

    .def("saveData", &gse::SEAppModule::saveData,
         R"eof(
Save results of state estimation calculation to data collection
objects
)eof")

    .def("hasConverged", &gse::SEAppModule::hasConverged)

    .def("preCheckMeasurements", &gse::SEAppModule::preCheckMeasurements,
         R"eof(
Perform pre-check for suspicious measurements
)eof")

    .def("debugMapper", &gse::SEAppModule::debugMapper,
         R"eof(
Perform a targeted check for bad measurements on Bus 8
)eof")

    .def("addVoltageLimitMeasurements",
         [](gse::SEAppModule& self, double vmin, double vmax, double deviation) {
           self.addVoltageLimitMeasurements(vmin, vmax, deviation);
         },
         py::arg("vmin") = 0.9,
         py::arg("vmax") = 1.1,
         py::arg("deviation") = 0.001,
         R"eof(
Add virtual measurements to enforce voltage magnitude constraints

Parameters:
     vmin (float): minimum allowed voltage magnitude
     vmax (float): maximum allowed voltage magnitude
     deviation (float): standard deviation for virtual measurements
)eof")

    .def("identifyPVBusConstraints",
         &gse::SEAppModule::identifyPVBusConstraints,
         R"eof(
Identify PV buses in the network and their connections Used to
properly handle voltage constraints at generator buses
)eof")

    .def("checkMeasurementConsistency",
         &gse::SEAppModule::checkMeasurementConsistency,
         R"eof( 
Check for potential measurement inconsistencies Identifies cases where
measurements may conflict with physical constraints
)eof")

    .def("handlePVBusVoltages", &gse::SEAppModule::handlePVBusVoltages,
         R"eof( 
Apply proper treatment for PV bus voltage measurements Ensures PV bus
voltages are treated as constraints rather than regular measurements
)eof")

    .def("handleVAMeasurements", &gse::SEAppModule::handleVAMeasurements,
         R"eof(
Apply special handling for voltage angle (VA) measurements Ensures
angle measurements at slack buses are treated as constraints
)eof")

    // Declared, but not implemented
    // .def("debugPrintMeasurements", &gse::SEAppModule::debugPrintMeasurements)

    .def("reportJacobianPerformance", &gse::SEAppModule::reportJacobianPerformance,
         R"eof(
Report Jacobian optimization performance statistics
)eof")

    ;
}


    
