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
// Added CA pybind11 bindings
// Updated March, 2026 by Yousu Chen 
// Added structured result bindings
// Updated September, 2026 by Yousu Chen
// -------------------------------------------------------------

#include "common.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"

namespace gpf = gridpack::powerflow;

// Compact float formatting for __repr__ strings. std::to_string always
// gives six decimal places, which prints a 1e-10 tolerance as "0.000000";
// %g keeps small and large magnitudes readable.
static std::string
repr_num(double v, int sig = 6)
{
  char buf[32];
  snprintf(buf, sizeof(buf), "%.*g", sig, v);
  return std::string(buf);
}

// -------------------------------------------------------------
// init_gridpack_pf
// Interface to powerflow application module
// -------------------------------------------------------------
void
init_gridpack_pf(py::module& gpm)
{
  py::module pfm =
    gpm.def_submodule("powerflow", "GridPACK power flow module");

  // Contingency type enum
  pfm.attr("Branch") = py::int_(static_cast<int>(gpf::Branch));
  pfm.attr("Generator") = py::int_(static_cast<int>(gpf::Generator));

  // Contingency struct for contingency analysis
  py::class_<gpf::Contingency>(pfm, "Contingency",
      "Container for contingency definition (line or generator)")
    .def(py::init<>())
    .def_readwrite("p_type", &gpf::Contingency::p_type,
        "Contingency type: 0=Generator, 1=Branch")
    .def_readwrite("p_name", &gpf::Contingency::p_name,
        "Contingency name")
    .def_readwrite("p_from", &gpf::Contingency::p_from,
        "From bus IDs for line contingencies")
    .def_readwrite("p_to", &gpf::Contingency::p_to,
        "To bus IDs for line contingencies")
    .def_readwrite("p_ckt", &gpf::Contingency::p_ckt,
        "Circuit IDs for line contingencies")
    // std::vector<bool> cannot use def_readwrite in pybind11 (proxy issue),
    // so use def_property with conversion to/from std::vector<int>.
    .def_property("p_saveLineStatus",
        [](const gpf::Contingency& self) {
          std::vector<int> v(self.p_saveLineStatus.begin(), self.p_saveLineStatus.end());
          return v;
        },
        [](gpf::Contingency& self, const std::vector<int>& v) {
          self.p_saveLineStatus.assign(v.begin(), v.end());
        },
        "Saved line status before contingency (as list of int: 1=true, 0=false)")
    .def_property("p_saveGenStatus",
        [](const gpf::Contingency& self) {
          std::vector<int> v(self.p_saveGenStatus.begin(), self.p_saveGenStatus.end());
          return v;
        },
        [](gpf::Contingency& self, const std::vector<int>& v) {
          self.p_saveGenStatus.assign(v.begin(), v.end());
        },
        "Saved generator status before contingency (as list of int: 1=true, 0=false)")
    .def_readwrite("p_busid", &gpf::Contingency::p_busid,
        "Bus IDs for generator contingencies")
    .def_readwrite("p_genid", &gpf::Contingency::p_genid,
        "Generator IDs for generator contingencies")
    ;

  // -----------------------------------------------------------------
  // Structured result types (gridpack/utilities/results_exporter.hpp),
  // populated by PFAppModule::solve / collectResults. Output-only: the
  // fields are read-only and there is no Python constructor.
  // -----------------------------------------------------------------

  py::class_<gpu::MismatchInfo>(pfm, "MismatchInfo",
      "Largest P and Q mismatch at one Newton iteration. Bus numbers are\n"
      "original PSS/E numbers. Globally reduced, so all ranks agree.")
    .def_readonly("maxPBus", &gpu::MismatchInfo::maxPBus, "bus number")
    .def_readonly("maxPMismatch", &gpu::MismatchInfo::maxPMismatch, "MW")
    .def_readonly("maxQBus", &gpu::MismatchInfo::maxQBus, "bus number")
    .def_readonly("maxQMismatch", &gpu::MismatchInfo::maxQMismatch, "MVAr")
    .def("__repr__", [](const gpu::MismatchInfo& m) {
        return "<MismatchInfo dP=" + repr_num(m.maxPMismatch) + " MW @ bus "
             + std::to_string(m.maxPBus) + ", dQ=" + repr_num(m.maxQMismatch)
             + " MVAr @ bus " + std::to_string(m.maxQBus) + ">";
      })
    ;

  py::class_<gpu::ConvergenceSummary>(pfm, "ConvergenceSummary",
      "Outcome of the most recent solve() or nl_solve(). Globally reduced\n"
      "inside solve(), so every rank sees the same summary and `converged`\n"
      "matches what solve() returned. Reads as zero before the first solve.")
    .def_readonly("converged", &gpu::ConvergenceSummary::converged,
        "True if tolerance was reached")
    .def_readonly("iterations", &gpu::ConvergenceSummary::iterations,
        "Newton iterations taken")
    .def_readonly("finalTolerance", &gpu::ConvergenceSummary::finalTolerance,
        "tolerance on the final iteration")
    .def_readonly("finalMismatch", &gpu::ConvergenceSummary::finalMismatch,
        "MismatchInfo for the final iteration")
    .def_readonly("perIteration", &gpu::ConvergenceSummary::perIteration,
        "MismatchInfo per iteration; each access copies the list")
    .def("__repr__", [](const gpu::ConvergenceSummary& c) {
        return std::string("<ConvergenceSummary ")
             + (c.converged ? "converged" : "DIVERGED")
             + " iterations=" + std::to_string(c.iterations)
             + " tol=" + repr_num(c.finalTolerance) + ">";
      })
    ;

  py::class_<gpu::BusResult>(pfm, "BusResult", "Solved quantities at one bus.")
    .def_readonly("busId", &gpu::BusResult::busId, "original bus number")
    .def_readonly("type", &gpu::BusResult::type, "1=PQ, 2=PV, 3=slack")
    .def_readonly("area", &gpu::BusResult::area, "area number")
    .def_readonly("zone", &gpu::BusResult::zone, "zone number")
    .def_readonly("baseKV", &gpu::BusResult::baseKV, "kV")
    .def_readonly("voltage", &gpu::BusResult::voltage, "pu")
    .def_readonly("angle", &gpu::BusResult::angle, "degrees")
    .def_readonly("pInjection", &gpu::BusResult::pInjection, "pGen - pLoad, MW")
    .def_readonly("qInjection", &gpu::BusResult::qInjection, "qGen - qLoad, MVAr")
    .def_readonly("pLoad", &gpu::BusResult::pLoad, "in-service load, MW")
    .def_readonly("qLoad", &gpu::BusResult::qLoad, "in-service load, MVAr")
    .def_readonly("pGen", &gpu::BusResult::pGen, "MW")
    .def_readonly("qGen", &gpu::BusResult::qGen, "MVAr")
    .def_readonly("shuntMvar", &gpu::BusResult::shuntMvar,
        "shunt MVAr seen by the Y-bus, fixed and switched combined")
    .def("__repr__", [](const gpu::BusResult& b) {
        return "<BusResult bus=" + std::to_string(b.busId)
             + " V=" + repr_num(b.voltage)
             + " pu, angle=" + repr_num(b.angle) + " deg>";
      })
    ;

  py::class_<gpu::BranchResult>(pfm, "BranchResult",
      "Solved flow on one branch circuit. PSS/E injection convention:\n"
      "pFrom and pTo both flow into the branch, so losses are their sum.")
    .def_readonly("fromBus", &gpu::BranchResult::fromBus, "original bus number")
    .def_readonly("toBus", &gpu::BranchResult::toBus, "original bus number")
    .def_readonly("circuitId", &gpu::BranchResult::circuitId, "circuit id")
    .def_readonly("pFrom", &gpu::BranchResult::pFrom, "MW at the from-bus")
    .def_readonly("qFrom", &gpu::BranchResult::qFrom, "MVAr at the from-bus")
    .def_readonly("pTo", &gpu::BranchResult::pTo, "MW at the to-bus")
    .def_readonly("qTo", &gpu::BranchResult::qTo, "MVAr at the to-bus")
    .def_readonly("pLoss", &gpu::BranchResult::pLoss, "pFrom + pTo, MW")
    .def_readonly("qLoss", &gpu::BranchResult::qLoss, "qFrom + qTo, MVAr")
    .def_readonly("mvaFrom", &gpu::BranchResult::mvaFrom, "MVA")
    .def_readonly("mvaTo", &gpu::BranchResult::mvaTo, "MVA")
    .def_readonly("rateA", &gpu::BranchResult::rateA,
        "rating A, MVA; zero when the source data gives no rating")
    .def_readonly("loadingPercent", &gpu::BranchResult::loadingPercent,
        "max(mvaFrom, mvaTo) / rateA * 100, or zero when rateA is zero")
    .def("__repr__", [](const gpu::BranchResult& b) {
        return "<BranchResult " + std::to_string(b.fromBus) + "-"
             + std::to_string(b.toBus) + " ckt " + b.circuitId
             + " P=" + repr_num(b.pFrom) + " MW, loading="
             + repr_num(b.loadingPercent) + "%>";
      })
    ;

  py::class_<gpu::GeneratorResult>(pfm, "GeneratorResult",
      "Solved dispatch of one generator.")
    .def_readonly("busId", &gpu::GeneratorResult::busId, "original bus number")
    .def_readonly("genId", &gpu::GeneratorResult::genId, "generator id")
    .def_readonly("pGen", &gpu::GeneratorResult::pGen, "MW")
    .def_readonly("qGen", &gpu::GeneratorResult::qGen, "MVAr")
    .def_readonly("qMax", &gpu::GeneratorResult::qMax, "upper Q limit, MVAr")
    .def_readonly("qMin", &gpu::GeneratorResult::qMin, "lower Q limit, MVAr")
    .def_readonly("voltageSetpoint", &gpu::GeneratorResult::voltageSetpoint,
        "scheduled voltage, pu")
    .def_readonly("status", &gpu::GeneratorResult::status, "1 in service, else 0")
    .def("__repr__", [](const gpu::GeneratorResult& g) {
        return "<GeneratorResult bus=" + std::to_string(g.busId)
             + " id=" + g.genId + " P=" + repr_num(g.pGen)
             + " MW, Q=" + repr_num(g.qGen) + " MVAr>";
      })
    ;

  py::class_<gpu::PowerFlowResults>(pfm, "PowerFlowResults",
      "Power flow solution from collectResults().\n\n"
      "Omitted from the lists: synthetic star buses from 3-winding\n"
      "transformers, their branches, and out-of-service circuits. Each\n"
      "list access copies, so bind it to a local name before looping.")
    .def_readonly("convergence", &gpu::PowerFlowResults::convergence,
        "ConvergenceSummary for the solve that produced these results")
    .def_readonly("buses", &gpu::PowerFlowResults::buses, "locally-owned buses")
    .def_readonly("branches", &gpu::PowerFlowResults::branches,
        "locally-owned, in-service branch circuits")
    .def_readonly("generators", &gpu::PowerFlowResults::generators,
        "generators on locally-owned buses")
    .def("__repr__", [](const gpu::PowerFlowResults& r) {
        return std::string("<PowerFlowResults ")
             + (r.convergence.converged ? "converged" : "DIVERGED")
             + " local: " + std::to_string(r.buses.size()) + " buses, "
             + std::to_string(r.branches.size()) + " branches, "
             + std::to_string(r.generators.size()) + " generators>";
      })
    ;

  py::enum_<gpu::ResultsExporter::Format>(pfm, "ExportFormat",
      "Output format for exportResults()")
    .value("JSON", gpu::ResultsExporter::JSON, "one '<basename>.json'")
    .value("CSV", gpu::ResultsExporter::CSV,
        "'<basename>_buses.csv' and _branches / _generators / _convergence")
    ;

  py::class_<gpf::PFAppModule> pfapp(pfm, "Powerflow");
  
  pfapp.doc() = ("Power flow application module");

  // minimal interface
  pfapp
    .def(py::init<>())
    .def("readNetwork",
         [](gpf::PFAppModule& self, gpu::Configuration& config, const int& idx,
            gpp::Communicator *comm) {

           // Default to the world communicator.  Contingency analysis
           // passes a task communicator instead so each task owns a
           // private copy of the network -- distributing tasks over ranks
           // that share one distributed network deadlocks in solve().
           gpp::Communicator net_comm;
           if (comm != NULL) net_comm = *comm;
           boost::shared_ptr<gpf::PFNetwork> pfnet( new gpf::PFNetwork(net_comm) );

           self.readNetwork(pfnet, &config, idx);
         }, py::arg("config"), py::arg("idx"), py::arg("comm") = nullptr,
         R"eof(
Read the network specified in the configuration.  

Parameters:
    config (gridpack.Configuration): power flow problem configuration usually read from an XML file
    idx (int): The index of the power flow problem in the configuration
    comm (gridpack.Communicator, optional): communicator to build the network
        on.  Defaults to the world communicator.  Pass a task communicator
        from :func:`Communicator.divide` to give each task its own copy.
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

#ifdef USE_OVERLOAD_CAST
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
#else
    .def("checkVoltageViolations",
         [](gpf::PFAppModule& self) -> py::object {
           bool flag(self.checkVoltageViolations());
           return py::cast(flag);
         },
         R"eof(
Check to see if there are any voltage violations in the network"

Returns:
    True if no violations found
)eof")
    .def("checkVoltageViolations",
         [](gpf::PFAppModule& self, const int& area) -> py::object {
           bool flag(self.checkVoltageViolations(area));
           return py::cast(flag);
         },
         R"eof(
Check to see if there are any voltage violations in an area

Parameters:
    area (int): area index

Returns:
    True if no violations found
)eof")
#endif
    .def("ignoreVoltageViolations", &gpf::PFAppModule::ignoreVoltageViolations,
         "Set \"ignore\" parameter on all buses with violations so that "
         "subsequent checks are not counted as violations")
    .def("clearVoltageViolations", &gpf::PFAppModule::clearVoltageViolations,
         "Clear \"ignore\" parameter on all buses")

#ifdef USE_OVERLOAD_CAST
    .def("checkLineOverloadViolations",
         py::overload_cast<>(&gpf::PFAppModule::checkLineOverloadViolations),
         "Check to see if there are any line overload violations in the network")
    .def("checkLineOverloadViolations",
         py::overload_cast<int>(&gpf::PFAppModule::checkLineOverloadViolations),
         "Check to see if there are any line overload violations in an area")
#else
    .def("checkLineOverloadViolations",
         [](gpf::PFAppModule& self) -> py::object {
           bool flag(self.checkLineOverloadViolations());
           return py::cast(flag);
         },
         R"eof(
Check to see if there are any line overload violations in the network.

Returns:
    True if no violations found
)eof")
    .def("checkLineOverloadViolations",
         [](gpf::PFAppModule& self, const int& area) -> py::object {
           bool flag(self.checkLineOverloadViolations(area));
           return py::cast(flag);
         },
         R"eof(
Check to see if there are any line overload violations in the network.

Parameters:
    area (int): area index

Returns:
    True if no violations found
)eof")
#endif
    .def("checkLineOverloadViolations",
         [](gpf::PFAppModule& self, std::vector<int> &bus1,
            std::vector<int> &bus2, std::vector<std::string> &tags) -> py::object {
           std::vector<bool> violations;
           bool flag =
             self.checkLineOverloadViolations(bus1, bus2, tags, violations);
           if (flag) {
             return py::cast(flag);
           } else {
             return py::cast(violations);
           }
         },
         R"eof(
Check for overload violations on specific lines

Parameters:
    bus1 : list of int
        original index of "from" bus for branch
    bus2 : list of int
        original index of "to" bus for branch
    tags : list of int
        line IDs for individual lines
  violations true if violation detected on branch, false otherwise


Returns:
    True if no violations found, list of flags for each line otherwise
)eof")
    
    // Declared but not implemented
    // .def("ignoreLineOverloadViolations",
    //      &gpf::PFAppModule::ignoreLineOverloadViolations,
    //      "Set \"ignore\" paramter on all lines with violations so "
    //      "that subsequent checks are not counted as violations")
    // .def("clearLineOverloadViolations",
    //      &gpf::PFAppModule::clearLineOverloadViolations,
    //      "Clear \"ignore\" parameter on all lines")

#ifdef USE_OVERLOAD_CAST
    .def("checkQlimViolations",
         py::overload_cast<>(&gpf::PFAppModule::checkQlimViolations),
         "Check to see if there are any Q limit violations in the network")
    .def("checkQlimViolations",
         py::overload_cast<int>(&gpf::PFAppModule::checkQlimViolations),
         "Check to see if there are any Q limit violations in an area")
#else
    .def("checkQlimViolations",
         [](gpf::PFAppModule& self) -> py::object {
           bool flag(self.checkQlimViolations());
           return py::cast(flag);
         },
         "Check to see if there are any Q limit violations in the network")
    .def("checkQlimViolations",
         [](gpf::PFAppModule& self, const int& area) -> py::object {
           bool flag(self.checkQlimViolations());
           return py::cast(flag);
         },
         "Check to see if there are any Q limit violations in an area")
#endif
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

    // Contingency analysis methods
    .def("setContingency", &gpf::PFAppModule::setContingency,
         R"eof(
Apply a contingency to the network (trip lines or generators).

Parameters:
    event (Contingency): contingency definition

Returns:
    True if contingency elements were found in network
)eof")
    .def("unSetContingency", &gpf::PFAppModule::unSetContingency,
         R"eof(
Restore network to pre-contingency state.

Parameters:
    event (Contingency): contingency to undo

Returns:
    True if restoration was successful
)eof")
    .def("writeCABus", &gpf::PFAppModule::writeCABus,
         "Write bus results for contingency analysis output")
    .def("writeCABranch", &gpf::PFAppModule::writeCABranch,
         "Write branch results for contingency analysis output")
    .def("checkSlackCapacity", &gpf::PFAppModule::checkSlackCapacity,
         "Check if slack bus generator is within capacity limits")
    .def("getIslandCount", &gpf::PFAppModule::getIslandCount,
         "Get number of islands in the network after contingency")
    .def("hasLoneBus", &gpf::PFAppModule::hasLoneBus,
         "Check if network has any isolated (lone) buses")

    ;

  // These query and modify the network's data collection objects
  pfapp
#ifdef USE_OVERLOAD_CAST
    .def("modifyDataCollectionGenParam",
         py::overload_cast<int, std::string, std::string, double>
         (&gpf::PFAppModule::modifyDataCollectionGenParam),
          "Modify (real) generator parameters in data collection for specified bus")
    .def("modifyDataCollectionGenParam",
         py::overload_cast<int, std::string, std::string, int>
         (&gpf::PFAppModule::modifyDataCollectionGenParam),
          "Modify (integer) generator parameters in data collection for specified bus")
#else
    .def("modifyDataCollectionGenParam",
         [](gpf::PFAppModule& self, int bus_id,
            std::string gen_id, std::string genParam,
            double value) -> py::object {
           bool flag =
             self.modifyDataCollectionGenParam(bus_id, gen_id, genParam, value);
           return py::cast(flag);
         }, 
         "Modify (real) generator parameters in data collection for specified bus")
    .def("modifyDataCollectionGenParam",
         [](gpf::PFAppModule& self, int bus_id,
            std::string gen_id, std::string genParam,
            int value) -> py::object {
           bool flag =
             self.modifyDataCollectionGenParam(bus_id, gen_id, genParam, value);
           return py::cast(flag);
         }, 
         "Modify (int) generator parameters in data collection for specified bus")
#endif

#ifdef USE_OVERLOAD_CAST
    .def("modifyDataCollectionLoadParam",
         py::overload_cast<int, std::string, std::string, double>
         (&gpf::PFAppModule::modifyDataCollectionLoadParam),
          "Modify (real) load parameters in data collection for specified bus")
    .def("modifyDataCollectionLoadParam",
         py::overload_cast<int, std::string, std::string, int>
         (&gpf::PFAppModule::modifyDataCollectionLoadParam),
          "Modify (integer) load parameters in data collection for specified bus")
#else
    .def("modifyDataCollectionLoadParam",
         [](gpf::PFAppModule& self, int bus_id,
            std::string gen_id, std::string genParam,
            double value) -> py::object {
           bool flag =
             self.modifyDataCollectionLoadParam(bus_id, gen_id, genParam, value);
           return py::cast(flag);
         },
         "Modify (real) load parameters in data collection for specified bus")
    .def("modifyDataCollectionLoadParam",
         [](gpf::PFAppModule& self, int bus_id,
            std::string gen_id, std::string genParam,
            int value) -> py::object {
           bool flag =
             self.modifyDataCollectionLoadParam(bus_id, gen_id, genParam, value);
           return py::cast(flag);
         }, 
         "Modify (integer) load parameters in data collection for specified bus")

#endif

#ifdef USE_OVERLOAD_CAST
    .def("modifyDataCollectionBusParam",
         py::overload_cast<int, std::string, double>
         (&gpf::PFAppModule::modifyDataCollectionBusParam),
          "Modify (real) parameters in data collection for specified bus")
    .def("modifyDataCollectionBusParam",
         py::overload_cast<int, std::string, int>
         (&gpf::PFAppModule::modifyDataCollectionBusParam),
          "Modify (integer) parameters in data collection for specified bus")
#else
    .def("modifyDataCollectionBusParam",
         [](gpf::PFAppModule& self, int bus_id, std::string busParam,
            double value) -> py::object {
           bool flag =
             self.modifyDataCollectionBusParam(bus_id, busParam, value);
           return py::cast(flag);
         },
         "Modify (real) parameters in data collection for specified bus")
    .def("modifyDataCollectionBusParam",
         [](gpf::PFAppModule& self, int bus_id, std::string busParam,
            int value) -> py::object {
           bool flag =
             self.modifyDataCollectionBusParam(bus_id, busParam, value);
           return py::cast(flag);
         }, 
         "Modify (integer) parameters in data collection for specified bus")
#endif

#ifdef USE_OVERLOAD_CAST
    .def("modifyDataCollectionBranchParam",
         py::overload_cast<int, int, std::string, std::string, double>
         (&gpf::PFAppModule::modifyDataCollectionBranchParam),
          "Modify (real) parameters in data collection for specified branch")
    .def("modifyDataCollectionBranchParam",
         py::overload_cast<int, int, std::string, std::string, int>
         (&gpf::PFAppModule::modifyDataCollectionBranchParam),
          "Modify (integer) parameters in data collection for specified branch")
#else
    .def("modifyDataCollectionBranchParam",
         [](gpf::PFAppModule& self, int bus1, int bus2, std::string ckt,
            std::string branchParam, double value) -> py::object {
           bool flag =
             self.modifyDataCollectionBranchParam(bus1, bus2, ckt,
                                                  branchParam, value);
           return py::cast(flag);
         },
         "Modify (real) parameters in data collection for specified branch")
    .def("modifyDataCollectionBranchParam",
         [](gpf::PFAppModule& self, int bus1, int bus2, std::string ckt,
            std::string branchParam, int value) -> py::object {
           bool flag =
             self.modifyDataCollectionBranchParam(bus1, bus2, ckt,
                                                  branchParam, value);
           return py::cast(flag);
         },
         "Modify (int) parameters in data collection for specified branch")
#endif
    
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

    .def("getConvergence", &gpf::PFAppModule::getConvergence,
         "ConvergenceSummary for the most recent solve(). nl_solve() does not\n"
       "populate it -- there it reads converged=false, iterations=0.")

    .def("collectResults", &gpf::PFAppModule::collectResults,
         R"eof(
Structured solution: buses, branches, generators, and convergence.

Values are read from the live network, so call it after solve(); modifying
data and re-solving, then calling again, returns the new solution.
)eof")

    .def("exportResults",
         [](gpf::PFAppModule& self, const std::string& basename,
            gpu::ResultsExporter::Format format) {
           self.exportResults(basename, format);
         },
         py::arg("basename"),
         py::arg("format") = gpu::ResultsExporter::JSON,
         R"eof(
Write the structured solution to JSON or CSV. `basename` is a stem, so
"run1.json" yields run1.json.json.

Only rank 0 writes, and only its own share, so on several ranks the file
is incomplete. Use it serially, or gather collectResults() yourself.
)eof")
    
    ;
    
}
