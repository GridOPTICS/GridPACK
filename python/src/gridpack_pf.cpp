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
         }, R"eof(
Read the network specified in the configuration.  

Parameters:
    config (gridpack.Configuration): power flow problem configuration usually read from an XML file
    idx (int): The index of the power flow problem in the configuration
)eof")
    .def("initialize", &gpf::PFAppModule::initialize,
         R"eof(
Initialize the power network. :func:`~gridpack.Powerflow.readNetwork` must be 
called before this.
)eof")
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
    .def("getPFSolutionSingleBus",
         [](gpf::PFAppModule& self, const int& bus_number) -> py::object {
           double bus_mag, bus_angle;
           bool flag(self.getPFSolutionSingleBus(bus_number, bus_mag, bus_angle));
           if (flag) {
             return py::make_tuple(bus_mag, bus_angle);
           } else {
             return py::cast<py::none>(Py_None);
           }
         },
         "get the power flow solution for the specific bus, vmag and v angle")
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
    .def("writeBus",
         [](gpf::PFAppModule& self, const std::string& signal) {
           self.writeBus(signal.c_str());
         },
         "Write bus results" )
    .def("writeBranch",
         [](gpf::PFAppModule& self, const std::string& signal) {
           self.writeBranch(signal.c_str());
         },
         "Write branch results" )
    .def("writeHeader",
         [](gpf::PFAppModule& self, const std::string& msg) {
           self.writeBranch(msg.c_str());
         },
         "Write the specified header string to output")
    .def("print",
         [](gpf::PFAppModule& self, const std::string& buf) {
           self.print(buf.c_str());
         },
         "Print an arbitrary string to output" )
    ;

  // Not implemented in Python, 
  // "writeBusString", "writeBranchString"

  pfapp
    .def("setVoltageLimits", &gpf::PFAppModule::setVoltageLimits,
         R"eof(
"Set voltage limits on all buses"

Parameters:
    Vmin (float): lower bound on voltages
    Vmax (float): upper bound on voltages
)eof")
         
    .def("checkVoltageViolations",
         py::overload_cast<>(&gpf::PFAppModule::checkVoltageViolations),
         R"eof(
Check to see if there are any voltage violations in the network"

Returns:
    True if no violations found
)eof")

    .def("checkVoltageViolations",
         py::overload_cast<int>(&gpf::PFAppModule::checkVoltageViolations),
         R"eof(
Check to see if there are any voltage violations in an area

Parameters:
    area (int): area index

Returns:
    True if no violations found
)eof")
    .def("ignoreVoltageViolations", &gpf::PFAppModule::ignoreVoltageViolations,
         "Set \"ignore\" parameter on all buses with violations so that "
         "subsequent checks are not counted as violations")
    .def("clearVoltageViolations", &gpf::PFAppModule::clearVoltageViolations,
         "Clear \"ignore\" parameter on all buses")

    .def("checkLineOverloadViolations",
         py::overload_cast<>(&gpf::PFAppModule::checkLineOverloadViolations),
         "Check to see if there are any line overload violations in the network")
    .def("checkLineOverloadViolations",
         py::overload_cast<int>(&gpf::PFAppModule::checkLineOverloadViolations),
         "Check to see if there are any line overload violations in an area")
    .def("checkLineOverloadViolations",
         py::overload_cast <
            std::vector<int>&, std::vector<int>&,
            std::vector<std::string>&, std::vector<bool>&
         >(&gpf::PFAppModule::checkLineOverloadViolations),
         "Check for overload violations on specific lines")

    // Declared but not implemented
    // .def("ignoreLineOverloadViolations",
    //      &gpf::PFAppModule::ignoreLineOverloadViolations,
    //      "Set \"ignore\" paramter on all lines with violations so "
    //      "that subsequent checks are not counted as violations")
    // .def("clearLineOverloadViolations",
    //      &gpf::PFAppModule::clearLineOverloadViolations,
    //      "Clear \"ignore\" parameter on all lines")

    
    .def("checkQlimViolations",
         py::overload_cast<>(&gpf::PFAppModule::checkQlimViolations),
         "Check to see if there are any Q limit violations in the network")
    .def("checkQlimViolations",
         py::overload_cast<int>(&gpf::PFAppModule::checkQlimViolations),
         "Check to see if there are any Q limit violations in an area")
    .def("clearQlimViolations", &gpf::PFAppModule::clearQlimViolations)
    .def("resetVoltages", &gpf::PFAppModule::resetVoltages)
    .def("scaleGeneratorRealPower", &gpf::PFAppModule::scaleGeneratorRealPower,
         R"eof(
Scale load power. 

Parameters:
    scale (float): factor to scale real power generation
    area (int): area index of area for scaling generation
    zone (int): zone index of zone for scaling generation

)eof")
    .def("scaleLoadPower", &gpf::PFAppModule::scaleLoadPower)
    .def("getTotalLoadRealPower", &gpf::PFAppModule::getTotalLoadRealPower)
    .def("getGeneratorMargins",
         [](gpf::PFAppModule& self, const int& area, const int& zone) -> py::object {
           double total(0.0), pmin(0.0), pmax(0.0);
           self.getGeneratorMargins(area, zone, &total, &pmin, &pmax);
           return py::make_tuple(total, pmin, pmax);
         },
         "Get the current real power generation and the maximum and minimum "
         "total power generation for all generators in a zone. If zone is less "
         "than 1 then return values for all generators in the area.")
    .def("resetPower", &gpf::PFAppModule::resetPower,
         "Reset power of loads and generators to original values")
    .def("writeRTPRDiagnostics",
         [](gpf::PFAppModule& self, const int& src_area, const int& src_zone,
            const int& load_area, const int& load_zone, const int& gen_scale,
            const int& load_scale, const std::string& file) {
           self.writeRTPRDiagnostics(src_area, src_zone,
                                     load_area, load_zone,
                                     gen_scale,load_scale, file.c_str());
         },
         "Write real time path rating diagnostics")
    .def("useRateB", &gpf::PFAppModule::useRateB,
         "Use rate B parameter for line overload violations")
         
    ;

  // These query and modify the network's data collection objects
  pfapp
    .def("modifyDataCollectionGenParam",
         py::overload_cast<int, std::string, std::string, double>
         (&gpf::PFAppModule::modifyDataCollectionGenParam),
          "Modify (real) generator parameters in data collection for specified bus")
    .def("modifyDataCollectionGenParam",
         py::overload_cast<int, std::string, std::string, int>
         (&gpf::PFAppModule::modifyDataCollectionGenParam),
          "Modify (integer) generator parameters in data collection for specified bus")

    .def("modifyDataCollectionLoadParam",
         py::overload_cast<int, std::string, std::string, double>
         (&gpf::PFAppModule::modifyDataCollectionLoadParam),
          "Modify (real) load parameters in data collection for specified bus")
    .def("modifyDataCollectionLoadParam",
         py::overload_cast<int, std::string, std::string, int>
         (&gpf::PFAppModule::modifyDataCollectionLoadParam),
          "Modify (integer) load parameters in data collection for specified bus")
    
    .def("modifyDataCollectionBusParam",
         py::overload_cast<int, std::string, double>
         (&gpf::PFAppModule::modifyDataCollectionBusParam),
          "Modify (real) parameters in data collection for specified bus")
    .def("modifyDataCollectionBusParam",
         py::overload_cast<int, std::string, int>
         (&gpf::PFAppModule::modifyDataCollectionBusParam),
          "Modify (integer) parameters in data collection for specified bus")

    .def("modifyDataCollectionBranchParam",
         py::overload_cast<int, int, std::string, std::string, double>
         (&gpf::PFAppModule::modifyDataCollectionBranchParam),
          "Modify (real) parameters in data collection for specified branch")
    .def("modifyDataCollectionBranchParam",
         py::overload_cast<int, int, std::string, std::string, int>
         (&gpf::PFAppModule::modifyDataCollectionBranchParam),
          "Modify (integer) parameters in data collection for specified branch")

    .def("getDataCollectionGenParamReal",
         [](gpf::PFAppModule& self, int bidx, std::string gen_id,
            std::string genParam) -> py::object {
           double result;
           bool ok = self.getDataCollectionGenParam(bidx, gen_id, genParam, &result);
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         },
         "Get (real) generator parameters in data collection for specified bus")

    .def("getDataCollectionGenParamInt",
         [](gpf::PFAppModule& self, int bidx, std::string gen_id,
            std::string genParam) -> py::object {
           int result;
           bool ok = self.getDataCollectionGenParam(bidx, gen_id, genParam, &result);
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         },
         "Get (integer) generator parameters in data collection for specified bus")
    
    .def("getDataCollectionLoadParamReal",
         [](gpf::PFAppModule& self, int bidx, std::string gen_id,
            std::string genParam) -> py::object {
           double result;
           bool ok = self.getDataCollectionLoadParam(bidx, gen_id, genParam, &result);
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         },
         "Get (real) load parameters in data collection for specified bus")

    .def("getDataCollectionLoadParamInt",
         [](gpf::PFAppModule& self, int bidx, std::string gen_id,
            std::string genParam) -> py::object {
           int result;
           bool ok = self.getDataCollectionLoadParam(bidx, gen_id, genParam, &result);
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         },
         "Get (integer) load parameters in data collection for specified bus")
    
    .def("getDataCollectionBusParamReal",
         [](gpf::PFAppModule& self, int bidx, std::string genParam) -> py::object {
           double result;
           bool ok = self.getDataCollectionBusParam(bidx, genParam, &result);
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         },
         "Get (real) load parameters in data collection for specified bus")

    .def("getDataCollectionBusParamInt",
         [](gpf::PFAppModule& self, int bidx, std::string genParam) -> py::object {
           int result;
           bool ok = self.getDataCollectionBusParam(bidx, genParam, &result);
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         },
         "Get (integer) load parameters in data collection for specified bus")
    
    .def("getDataCollectionBranchParamReal",
         [](gpf::PFAppModule& self, int bidx1, int bidx2, std::string ckt, 
            std::string genParam) -> py::object {
           double result;
           bool ok =
             self.getDataCollectionBranchParam(bidx1, bidx2, ckt, genParam, &result);
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         },
         "Get (real) load parameters in data collection for specified bus")

    .def("getDataCollectionBranchParamInt",
         [](gpf::PFAppModule& self, int bidx1, int bidx2, std::string ckt, 
            std::string genParam) -> py::object {
           int result;
           bool ok =
             self.getDataCollectionBranchParam(bidx1, bidx2, ckt, genParam, &result);
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         },
         "Get (integer) load parameters in data collection for specified bus")
    
    ;
    
}
