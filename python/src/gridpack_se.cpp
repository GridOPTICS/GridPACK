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
// Last Change: 2025-09-22 14:59:20 d3g096
// -------------------------------------------------------------

#include <pybind11/stl_bind.h>
#include <cstring>
#include <stdexcept>
#include "common.hpp"
#include "gridpack/applications/modules/state_estimation/se_app_module.hpp"

namespace gse = gridpack::state_estimation;

PYBIND11_MAKE_OPAQUE(std::vector<gse::Measurement>)

// -------------------------------------------------------------
// init_gridpack_se
// Interface to the GridPACK state estimation application module
// -------------------------------------------------------------
void
init_gridpack_se(py::module& gpm)
{
  py::module sem =
    gpm.def_submodule("state_estimation", "GridPACK state estimation module");


  // -------------------------------------------------------------
  // gridpack.state_estimation.Measurement
  // -------------------------------------------------------------
  py::class_<gse::Measurement> meas(sem, "Measurement");

  meas.doc() = ("Container for state estimation measurement data");

  meas
    .def(py::init<>())
    .def_readwrite("busid", &gse::Measurement::p_busid)
    .def_readwrite("fbusid", &gse::Measurement::p_fbusid)
    .def_readwrite("tbusid", &gse::Measurement::p_tbusid)
    .def_readwrite("value", &gse::Measurement::p_value)
    .def_readwrite("deviation", &gse::Measurement::p_deviation)
    .def("type", [](gse::Measurement &self) -> std::string {
      std::string result(self.p_type);
      return result;
    }, R"eof(
Get measurement type.
)eof")
    .def("type", [](gse::Measurement &self, const std::string& code) {
      if (code.length() > 4) {
        std::string msg("Invalid Measurement type: ");
        msg += code;
        throw std::runtime_error(msg);
      }
      strncpy(self.p_type, code.c_str(), 4);
    }, R"eof(
Set measurment type, one of "VA", "VM", "VUL", "VLL", "PIJ", "QIJ"
)eof")
    .def("ckt", [](gse::Measurement &self) -> std::string {
      std::string result(self.p_ckt);
      return result;
    }, R"eof(
Get measurement ckt.
)eof")
    .def("ckt", [](gse::Measurement &self, const std::string& code) {
      if (code.length() > 3) {
        std::string msg("Invalid Measurement ckt: ");
        msg += code;
        throw std::runtime_error(msg);
      }
      strncpy(self.p_ckt, code.c_str(), 3);
    }, R"eof(
Set measurment ckt"
)eof")
    ;


  // -------------------------------------------------------------
  // gridpack.state_estimation.MeasurementVector
  // -------------------------------------------------------------

  py::bind_vector< std::vector< gse::Measurement > >(sem, "MeasurementVector");

  // -------------------------------------------------------------
  // gridpack.state_estimation.SEApp
  // -------------------------------------------------------------
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

    .def("getMeasurements",
         [](gse::SEAppModule& self, gpu::Configuration& config) ->
         std::vector<gse::Measurement> {
           gpu::Configuration::CursorPtr cursor;
           cursor = config.getCursor("Measurements");
           gpu::Configuration::ChildCursors mcursors;
           cursor->children(mcursors);
           std::vector<gse::Measurement> measures;
           measures = self.getMeasurements(mcursors);
           return(measures);
         }, R"eof(
Get measurments from an open Configuration. 
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


    
