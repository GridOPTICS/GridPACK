// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: gridpack_ds.hpp
//
// Interface to Dynamic Simulation and HADREC applications
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

#include <gridpack/applications/modules/hadrec/hadrec_app_module.hpp>
#include <gridpack/applications/modules/powerflow/pf_app_module.hpp>

namespace gph = gridpack::hadrec;
namespace gpds = gridpack::dynamic_simulation;

// Some pybind11 magic for a vector of Event
PYBIND11_MAKE_OPAQUE( std::vector< gpds::Event > )

// -------------------------------------------------------------
// init_gridpack_ds
// Initialize Dynamic Simulation application inteface
// -------------------------------------------------------------
void
init_gridpack_ds(py::module& gpm)
{
  // -------------------------------------------------------------
  // gridpack.dynamic_simulation module
  // -------------------------------------------------------------
  py::module dsm =
    gpm.def_submodule("dynamic_simulation",
                      "GridPACK Dynamic Simulation Application module");

  // -------------------------------------------------------------
  // gridpack.dynamic_simulation.EventVector
  // -------------------------------------------------------------
  py::bind_vector< std::vector< gpds::Event > >(dsm, "EventVector");

  // -------------------------------------------------------------
  // gridpack.dynamic_simulation.Event
  // -------------------------------------------------------------
  py::class_<gpds::Event>(dsm, "Event")
    .def(py::init<>())
    .def_readwrite("start", &gpds::Event::start)
    .def_readwrite("end", &gpds::Event::end)
    .def_readwrite("step", &gpds::Event::step)
    .def_readwrite("tag", &gpds::Event::tag)
    .def_readwrite("isGenerator", &gpds::Event::isGenerator)
    .def_readwrite("isBus", &gpds::Event::isBus)
    .def_readwrite("bus_idx", &gpds::Event::bus_idx)
    .def_readwrite("isLine", &gpds::Event::isLine)
    .def_readwrite("from_idx", &gpds::Event::from_idx)
    .def_readwrite("to_idx", &gpds::Event::to_idx)
    ;

  // -------------------------------------------------------------
  // gridpack.dynamic_simulation.DSFullApp
  // -------------------------------------------------------------
  py::class_<gpds::DSFullApp> dsapp(dsm, "DSFullApp");
  dsapp
    .def(py::init<>())
    .def("solvePowerFlowBeforeDynSimu",
         [](gpds::DSFullApp& self, const std::string& inputfile, const int& pf_idx) {
           self.solvePowerFlowBeforeDynSimu(inputfile.c_str(), pf_idx);
         })
    // Alternative init that uses existing Configuration singleton and
    // communicator, avoiding MPI_Comm_dup conflicts when Environment
    // is already initialized (e.g., after power flow or CA).
    .def("initFromConfig",
         [](gpds::DSFullApp& self, gpu::Configuration& config, const int& pf_idx) {
           gpp::Communicator world;

           // Run power flow
           boost::shared_ptr<gridpack::powerflow::PFNetwork>
             pf_network(new gridpack::powerflow::PFNetwork(world));
           gridpack::powerflow::PFAppModule pf_app;
           pf_app.readNetwork(pf_network, &config, pf_idx);
           pf_app.initialize();
           // Same solver selection and report as
           // DSFullApp::solvePowerFlowBeforeDynSimu.
           gpu::Configuration::CursorPtr cursor =
             config.getCursor("Configuration.Powerflow");
           bool useNonLinear = false;
           useNonLinear = cursor->get("UseNonLinear", useNonLinear);
           if (useNonLinear) {
             pf_app.nl_solve();
           } else {
             pf_app.solve();
           }
           pf_app.write();
           pf_app.saveData();

           // Create DS network and transfer PF results
           boost::shared_ptr<gpds::DSFullNetwork>
             ds_network(new gpds::DSFullNetwork(world));
           pf_network->clone<gpds::DSFullBus, gpds::DSFullBranch>(ds_network);
           self.transferPFtoDS(pf_network, ds_network);
           self.setNetwork(ds_network, &config);
         },
         R"eof(
Initialize dynamic simulation from an already-open Configuration.

Uses the existing MPI communicator instead of creating a new one,
avoiding MPI_Comm_dup conflicts when gridpack.Environment() is
already initialized.

Parameters:
    config (gridpack.Configuration): Already-open configuration
    pf_idx (int): Power flow configuration index (-1 for default)
)eof")
    .def("readGenerators", &gpds::DSFullApp::readGenerators,
         py::arg("ds_idx") = -1)
    .def("readSequenceData", &gpds::DSFullApp::readSequenceData)
    .def("initialize", &gpds::DSFullApp::initialize)
    .def("reload", &gpds::DSFullApp::reload)
    .def("reset", &gpds::DSFullApp::reset)
    .def("solve", &gpds::DSFullApp::solve)
    .def("solvePreInitialize", &gpds::DSFullApp::solvePreInitialize)
    .def("setup", &gpds::DSFullApp::setup)
    .def("executeOneSimuStep", &gpds::DSFullApp::executeOneSimuStep)
#ifdef USE_OVERLOAD_CAST
    .def("run", py::overload_cast<>(&gpds::DSFullApp::run))
    .def("run", py::overload_cast<double>(&gpds::DSFullApp::run))
#else
    .def("run", [](gpds::DSFullApp& self) {self.run();})
    .def("run", [](gpds::DSFullApp& self, double tend) {self.run(tend);})
#endif
    .def("scatterInjectionLoad", &gpds::DSFullApp::scatterInjectionLoad)
    .def("scatterInjectionLoadNew", &gpds::DSFullApp::scatterInjectionLoadNew)
    .def("scatterInjectionLoadNew_compensateY",
         &gpds::DSFullApp::scatterInjectionLoadNew_compensateY)
    .def("scatterInjectionLoadNew_Norton",
         &gpds::DSFullApp::scatterInjectionLoadNew_Norton)
    .def("scatterInjectionLoadNewConstCur",
         &gpds::DSFullApp::scatterInjectionLoadNewConstCur)
    .def("applyLoadShedding", &gpds::DSFullApp::applyLoadShedding)
    .def("applyConstYLoadShedding", &gpds::DSFullApp::applyConstYLoadShedding)
    .def("setWideAreaControlSignal", &gpds::DSFullApp::setWideAreaControlSignal)
    .def("applyGFIAdjustment", &gpds::DSFullApp::applyGFIAdjustment)
    .def("applyConstYLoad_Change_P", &gpds::DSFullApp::applyConstYLoad_Change_P)
    .def("clearConstYLoad_Change_P", &gpds::DSFullApp::clearConstYLoad_Change_P)
    .def("applyConstYLoad_Change_Q", &gpds::DSFullApp::applyConstYLoad_Change_Q)
    .def("clearConstYLoad_Change_Q", &gpds::DSFullApp::clearConstYLoad_Change_Q)
    .def("setConstYLoadtoZero_P", &gpds::DSFullApp::setConstYLoadtoZero_P)
    .def("setConstYLoadtoZero_Q", &gpds::DSFullApp::setConstYLoadtoZero_Q)
    .def("setConstYLoadImpedance", &gpds::DSFullApp::setConstYLoadImpedance)
    .def("applyGeneratorTripping", &gpds::DSFullApp::applyGeneratorTripping)
#ifdef USE_OVERLOAD_CAST
    .def("setLineTripAction",
         py::overload_cast<int, int, std::string>(&gpds::DSFullApp::setLineTripAction))
    .def("setLineTripAction",
         py::overload_cast<int>(&gpds::DSFullApp::setLineTripAction))
#else
#endif
    .def("clearLineTripAction", &gpds::DSFullApp::clearLineTripAction)
    .def("isDynSimuDone", &gpds::DSFullApp::isDynSimuDone)
    .def("write",
         [](gpds::DSFullApp& self, const std::string& signal) {
           self.write(signal.c_str());
         })
    ;

  dsapp
    .def("setObservations",
         [] (gpds::DSFullApp& self, ConfigurationCursorWrapper& cursor) {
           self.setObservations(cursor.the_cursor);
         })
    .def("getEvents",
         [](gpds::DSFullApp& self) {
           return self.getEvents();
         },
         py::return_value_policy::copy)
    .def("getEvents",
         [](gpds::DSFullApp& self, ConfigurationCursorWrapper c) {
           return (self.getEvents(c.the_cursor));
         },
         py::return_value_policy::copy)
    .def("setEvent", &gpds::DSFullApp::setEvent)
    ;

  dsapp
    .def("setGeneratorWatch",
         [](gpds::DSFullApp& self) { self.setGeneratorWatch(); })
    .def("setGeneratorWatch",
         [](gpds::DSFullApp& self, const std::string& filename) {
           self.setGeneratorWatch(filename.c_str());
         })
    .def("setGeneratorWatch",
         [](gpds::DSFullApp& self, std::vector<int>& buses,
            std::vector<std::string>& tags, bool writeFile) {
           self.setGeneratorWatch(buses, tags, writeFile);
         },
         py::arg("buses") = std::vector<int>(),
         py::arg("tags") = std::vector<std::string>(),
         py::arg("writeFile") = true
         )
    .def("setGeneratorWatch",
         [](gpds::DSFullApp& self, const std::string& filename,
            ConfigurationCursorWrapper& cursor) {
           self.setGeneratorWatch(filename.c_str(), cursor.the_cursor);
         })
    .def("setGeneratorWatch",
         [](gpds::DSFullApp& self, ConfigurationCursorWrapper& cursor) {
           self.setGeneratorWatch(cursor.the_cursor);
         })
    .def("setLoadWatch", &gpds::DSFullApp::setLoadWatch)
    ;

  dsapp
    .def("open",
         [](gpds::DSFullApp& self, const std::string& filename) {
           self.open(filename.c_str());
         })
    .def("close", &gpds::DSFullApp::close)
    .def("print",
         [](gpds::DSFullApp& self, const std::string& buf) {
           self.print(buf.c_str());
         })
    .def("isSecure", &gpds::DSFullApp::isSecure)
    .def("saveTimeSeries", &gpds::DSFullApp::saveTimeSeries)
    ;

  dsapp
    .def("getTimeSeriesMap",
         &gpds::DSFullApp::getTimeSeriesMap,
         py::return_value_policy::copy)
    .def("getGeneratorTimeSeries",
         &gpds::DSFullApp::getGeneratorTimeSeries,
         py::return_value_policy::copy)
    .def("getListWatchedGenerators",
         [](gpds::DSFullApp& self) -> py::object {
           std::vector<int> bus_ids;
           std::vector<std::string> gen_ids;
           self.getListWatchedGenerators(bus_ids, gen_ids);
           return py::make_tuple(bus_ids, gen_ids);
         })
    .def("frequencyOK", &gpds::DSFullApp::frequencyOK)
    .def("scaleGeneratorRealPower", &gpds::DSFullApp::scaleGeneratorRealPower)
    .def("scaleLoadPower", &gpds::DSFullApp::scaleLoadPower)
    .def("getTotalLoadRealPower", &gpds::DSFullApp::getTotalLoadRealPower)
    .def("getGeneratorMargins",
         [](gpds::DSFullApp& self, int area, int zone) -> py::object {
           double total, pmin, pmax;
           self.getGeneratorMargins(area, zone, &total, &pmin, &pmax);
           return py::make_tuple(total, pmin, pmax);
         })
    .def("resetPower", &gpds::DSFullApp::resetPower)
    .def("updateData", &gpds::DSFullApp::updateData)
    .def("writeRTPRDiagnostics",
         [](gpds::DSFullApp& self,
            int src_area, int src_zone, int load_area,
            int load_zone, double gen_scale, double load_scale,
            std::string file) {
           self.writeRTPRDiagnostics(src_area, src_zone, load_area,
                                     load_zone, gen_scale, load_scale,
                                     file.c_str());
         })
    .def("getFrequencyFailures", &gpds::DSFullApp::getFrequencyFailures,
         py::return_value_policy::copy)
    .def("setFrequencyMonitoring", &gpds::DSFullApp::setFrequencyMonitoring)
    ;

  dsapp
    .def("getObservationLists",
         [](gpds::DSFullApp& self) -> py::object {
           std::vector<int> obs_genBus;
           std::vector<std::string> obs_genIDs;
           std::vector<int> obs_loadBuses;
           std::vector<std::string> obs_loadIDs;
           std::vector<int> obs_busIDs;
           self.getObservationLists(obs_genBus, obs_genIDs,
                                    obs_loadBuses, obs_loadIDs, obs_busIDs);
           return py::make_tuple(obs_genBus, obs_genIDs,
                                 obs_loadBuses, obs_loadIDs, obs_busIDs);
         })
    .def("getObservationLists_withBusFreq",
         [](gpds::DSFullApp& self) -> py::object {
           std::vector<int> obs_genBus;
           std::vector<std::string> obs_genIDs;
           std::vector<int> obs_loadBuses;
           std::vector<std::string> obs_loadIDs;
           std::vector<int> obs_busIDs;
           std::vector<int> busfreqIDs;
           self.getObservationLists_withBusFreq(obs_genBus, obs_genIDs,
                                                obs_loadBuses, obs_loadIDs,
                                                obs_busIDs, busfreqIDs);
           return py::make_tuple(obs_genBus, obs_genIDs,
                                 obs_loadBuses, obs_loadIDs,
                                 obs_busIDs, busfreqIDs);
         })
    .def("getObservations",
         [](gpds::DSFullApp& self) -> py::object {
           std::vector<double> vMag, vAng, rSpd, rAng, genP, genQ, fOnline;
           self.getObservations(vMag, vAng, rSpd, rAng, genP, genQ, fOnline);
           return py::make_tuple(vMag, vAng, rSpd, rAng, genP, genQ, fOnline);
         })
    .def("getObservations_withBusFreq", 
         [](gpds::DSFullApp& self) -> py::object {
           std::vector<double> vMag, vAng, rSpd, rAng, genP, genQ, fOnline, busfreq;
           self.getObservations_withBusFreq(vMag, vAng, rSpd, rAng,
                                            genP, genQ, fOnline, busfreq);
           return py::make_tuple(vMag, vAng, rSpd, rAng,
                                 genP, genQ, fOnline, busfreq);
         })
    ;         

  dsapp
    .def("getBusTotalLoadPower",
         [](gpds::DSFullApp& self, const int& busid) -> py::object {
           double pg, qg;
           bool flag;
           flag = self.getBusTotalLoadPower(busid, pg, qg);
           if (flag) {
             return py::make_tuple(pg, qg);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getGeneratorPower",
         [](gpds::DSFullApp& self, int bus_id, std::string gen_id) -> py::object {
           double pg, qg;
           bool flag(self.getGeneratorPower(bus_id, gen_id, pg, qg));
           if (flag) {
             return py::make_tuple(pg, qg);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getZoneLoads",
         [](const gpds::DSFullApp& self) -> py::object {
           std::vector<double> load_p, load_q;
           std::vector<int> zone_id;
           self.getZoneLoads(load_p, load_q, zone_id);
           return py::make_tuple(load_p, load_q, zone_id);
         })
    .def("getZoneGeneratorPower",
         [](const gpds::DSFullApp& self) -> py::object {
           std::vector<double> generator_p, generator_q;
           std::vector<int> zone_id;
           self.getZoneGeneratorPower(generator_p, generator_q, zone_id);
           return py::make_tuple(generator_p, generator_q, zone_id);
         })
    ;

  dsapp
#ifdef USE_OVERLOAD_CAST
    .def("modifyDataCollectionGenParam",
         py::overload_cast<int, std::string, std::string, double>
         (&gpds::DSFullApp::modifyDataCollectionGenParam))
    .def("modifyDataCollectionGenParam",
         py::overload_cast<int, std::string, std::string, int>
         (&gpds::DSFullApp::modifyDataCollectionGenParam))
#else
    .def("modifyDataCollectionGenParam",
         [](gpds::DSFullApp& self, int bus_id,
            std::string gen_id, std::string genParam,
            double value) -> py::object {
           bool flag =
             self.modifyDataCollectionGenParam(bus_id, gen_id, genParam, value);
           return py::cast(flag);
         })
    .def("modifyDataCollectionGenParam",
         [](gpds::DSFullApp& self, int bus_id,
            std::string gen_id, std::string genParam,
            int value) -> py::object {
           bool flag =
             self.modifyDataCollectionGenParam(bus_id, gen_id, genParam, value);
           return py::cast(flag);
         })
#endif

#ifdef USE_OVERLOAD_CAST
    .def("modifyDataCollectionLoadParam",
         py::overload_cast<int, std::string, std::string, double>
         (&gpds::DSFullApp::modifyDataCollectionLoadParam))
    .def("modifyDataCollectionLoadParam",
         py::overload_cast<int, std::string, std::string, int>
         (&gpds::DSFullApp::modifyDataCollectionLoadParam))
#else
    .def("modifyDataCollectionLoadParam",
         [](gpds::DSFullApp& self, int bus_id, std::string load_id,
            std::string loadParam, double value) -> py::object {
           bool flag =
             self.modifyDataCollectionLoadParam(bus_id, load_id,
                                                loadParam, value);
           return py::cast(flag);
         })
    .def("modifyDataCollectionLoadParam",
         [](gpds::DSFullApp& self, int bus_id, std::string load_id,
            std::string loadParam, int value) -> py::object {
           bool flag =
             self.modifyDataCollectionLoadParam(bus_id, load_id,
                                                loadParam, value);
           return py::cast(flag);
         })
#endif

#ifdef USE_OVERLOAD_CAST
    .def("modifyDataCollectionBusParam",
         py::overload_cast<int, std::string, double>
         (&gpds::DSFullApp::modifyDataCollectionBusParam))
    .def("modifyDataCollectionBusParam",
         py::overload_cast<int, std::string, int>
         (&gpds::DSFullApp::modifyDataCollectionBusParam))
#else
    .def("modifyDataCollectionBusParam",
         [](gpds::DSFullApp& self, int bus_id,
            std::string busParam, double value) -> py::object {
           bool flag =
             self.modifyDataCollectionBusParam(bus_id, busParam, value);
           return py::cast(flag);
         })
    .def("modifyDataCollectionBusParam",
         [](gpds::DSFullApp& self, int bus_id,
            std::string busParam, int value) -> py::object {
           bool flag =
             self.modifyDataCollectionBusParam(bus_id, busParam, value);
           return py::cast(flag);
         })
#endif
    
    ;   

  dsapp
    .def("setState", &gpds::DSFullApp::setState)
    .def("getState",
         [](gpds::DSFullApp& self, const int& bus_id,
            const std::string& dev_id, const std::string& device,
            const std::string& name) -> py::object {
           bool status;
           double value;
           status = self.getState(bus_id, dev_id, device, name, &value);

           // Value is returned if successful.
           if (status) {
             return  py::cast(value);
           }

           // Otherwise, None is returned. 
           return py::object();
         })
    .def("getTimeStep", &gpds::DSFullApp::getTimeStep)
    .def("setTimeStep", &gpds::DSFullApp::setTimeStep)
    .def("setFinalTime", &gpds::DSFullApp::setFinalTime)
    .def("getFinalTime", &gpds::DSFullApp::getFinalTime)
    .def("getCurrentTime", &gpds::DSFullApp::getCurrentTime)
    ;

  dsapp
    .def("totalBuses", &gpds::DSFullApp::totalBuses)
    .def("totalBranches", &gpds::DSFullApp::totalBranches)
    .def("getConnectedBranches", &gpds::DSFullApp::getConnectedBranches,
         py::return_value_policy::copy)
    .def("getBranchEndpoints",
         [](gpds::DSFullApp& self, const int& oidx) {
           int fbus, tbus;
           self.getBranchEndpoints(oidx, &fbus, &tbus);
           return py::make_tuple(fbus, tbus);
         })
#ifdef USE_OVERLOAD_CAST
    .def("numGenerators",
         py::overload_cast<>(&gpds::DSFullApp::numGenerators, py::const_))
    .def("numGenerators",
         py::overload_cast<const int&>(&gpds::DSFullApp::numGenerators, py::const_))
#else
    .def("numGenerators",
         [](const gpds::DSFullApp& self) -> py::object {
           int n(self.numGenerators());
           return py::cast(n);
         })
    .def("numGenerators",
         [](const gpds::DSFullApp& self, const int& bus_idx) -> py::object {
           int n(self.numGenerators(bus_idx));
           return py::cast(n);
         })
#endif

#ifdef USE_OVERLOAD_CAST
    .def("numLoads",
         py::overload_cast<>(&gpds::DSFullApp::numLoads, py::const_))
    .def("numLoads",
         py::overload_cast<const int&>(&gpds::DSFullApp::numLoads, py::const_))
#else
    .def("numLoads",
         [](const gpds::DSFullApp& self) -> py::object {
           int n(self.numLoads());
           return py::cast(n);
         })
    .def("numLoads",
         [](const gpds::DSFullApp& self, const int& bus_idx) -> py::object {
           int n(self.numLoads(bus_idx));
           return py::cast(n);
         })
#endif

#ifdef USE_OVERLOAD_CAST    
    .def("numStorage",
         py::overload_cast<>(&gpds::DSFullApp::numStorage, py::const_))
    .def("numStorage",
         py::overload_cast<const int&>(&gpds::DSFullApp::numStorage, py::const_))
#else
    .def("numStorage",
         [](const gpds::DSFullApp& self) -> py::object {
           int n(self.numStorage());
           return py::cast(n);
         })
    .def("numStorage",
         [](const gpds::DSFullApp& self, const int& bus_idx) -> py::object {
           int n(self.numStorage(bus_idx));
           return py::cast(n);
         })
#endif

#ifdef USE_OVERLOAD_CAST    
    .def("numLines",
         py::overload_cast<>(&gpds::DSFullApp::numLines, py::const_))
    .def("numLines",
         py::overload_cast<const int&>(&gpds::DSFullApp::numLines, py::const_))
#else
    .def("numLines",
         [](const gpds::DSFullApp& self) -> py::object {
           int n(self.numLines());
           return py::cast(n);
         })
    .def("numLines",
         [](const gpds::DSFullApp& self, const int& branch_idx) -> py::object {
           int n(self.numLines(branch_idx));
           return py::cast(n);
         })

#endif
    
    .def("getBusInfoInt", 
         [](gpds::DSFullApp& self, const int& bus_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           int result;
           bool ok = self.getBusInfo(bus_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBusInfoBool",
         [](gpds::DSFullApp& self, const int& bus_idx,
            const std::string& name, const std::optional<int> dev_idx = -1) -> py::object {
           bool result;
           bool ok = self.getBusInfo(bus_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getBusInfoReal",
         [](gpds::DSFullApp& self, const int& bus_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           double result;
           bool ok = self.getBusInfo(bus_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBusInfoString", 
         [](gpds::DSFullApp& self, const int& bus_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           std::string result;
           bool ok = self.getBusInfo(bus_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBranchInfoInt", 
         [](gpds::DSFullApp& self, const int& branch_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           int result;
           bool ok = self.getBranchInfo(branch_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBranchInfoBool",
         [](gpds::DSFullApp& self, const int& branch_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           bool result;
           bool ok = self.getBranchInfo(branch_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBranchInfoReal",
         [](gpds::DSFullApp& self, const int& branch_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           double result;
           bool ok = self.getBranchInfo(branch_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBranchInfoString", 
         [](gpds::DSFullApp& self, const int& branch_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           std::string result;
           bool ok = self.getBranchInfo(branch_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    ;

  
}


// -------------------------------------------------------------
// init_gridpack_hadrec
// Initialize the HADREC application interface
// -------------------------------------------------------------
void
init_gridpack_hadrec(py::module & gpm)
{
  py::module hadm =
    gpm.def_submodule("hadrec", "GridPACK HADREC Application module");
  
  // -------------------------------------------------------------
  // gridpack.hadrec.Action
  // -------------------------------------------------------------
  py::class_<gph::HADRECAction>(hadm, "Action")
    .def(py::init<>())
    .def_readwrite("actiontype", &gph::HADRECAction::actiontype)
    .def_readwrite("bus_number", &gph::HADRECAction::bus_number)
    .def_readwrite("componentID", &gph::HADRECAction::componentID)
    .def_readwrite("percentage", &gph::HADRECAction::percentage)
    .def_readwrite("brch_from_bus_number", &gph::HADRECAction::brch_from_bus_number)
    .def_readwrite("brch_to_bus_number", &gph::HADRECAction::brch_to_bus_number)
    .def_readwrite("branch_ckt", &gph::HADRECAction::branch_ckt)
    ;

  // -------------------------------------------------------------
  // gridpack.hadrec.Module
  // -------------------------------------------------------------
  py::class_<gph::HADRECAppModule> hadapp(hadm, "Module");
  hadapp
    .def(py::init<>())
    .def("transferPFtoDS", &gph::HADRECAppModule::transferPFtoDS)
    .def("executeDynSimuOneStep", &gph::HADRECAppModule::executeDynSimuOneStep)
    .def("updateData", &gph::HADRECAppModule::updateData)
    .def("isDynSimuDone",  &gph::HADRECAppModule::isDynSimuDone)
    .def("applyAction", &gph::HADRECAppModule::applyAction)
    .def("getObservations", &gph::HADRECAppModule::getObservations,
         py::return_value_policy::copy)
    ;

  // These are some network topology/analytics query methods
  hadapp
    .def("totalBuses", &gph::HADRECAppModule::totalBuses)
    .def("totalBranches", &gph::HADRECAppModule::totalBranches)
    .def("getConnectedBranches", &gph::HADRECAppModule::getConnectedBranches,
         py::return_value_policy::copy)
    .def("getBranchEndpoints",
         [](gph::HADRECAppModule& self, const int& oidx) {
           int fbus, tbus;
           self.getBranchEndpoints(oidx, &fbus, &tbus);
           return py::make_tuple(fbus, tbus);
         })

#ifdef USE_OVERLOAD_CAST
    .def("numGenerators",
         py::overload_cast<>(&gph::HADRECAppModule::numGenerators, py::const_))
    .def("numGenerators",
         py::overload_cast<const int&>(&gph::HADRECAppModule::numGenerators, py::const_))
#else
    .def("numGenerators",
         [](const gph::HADRECAppModule& self) -> py::object {
           int n(self.numGenerators());
           return py::cast(n);
         })
    .def("numGenerators",
         [](const gph::HADRECAppModule& self, const int& bus_idx) -> py::object {
           int n(self.numGenerators(bus_idx));
           return py::cast(n);
         })
#endif

#ifdef USE_OVERLOAD_CAST
    .def("numLoads",
         py::overload_cast<>(&gph::HADRECAppModule::numLoads, py::const_))
    .def("numLoads",
         py::overload_cast<const int&>(&gph::HADRECAppModule::numLoads, py::const_))
#else
    .def("numLoads",
         [](const gph::HADRECAppModule& self) -> py::object {
           int n(self.numLoads());
           return py::cast(n);
         })
    .def("numLoads",
         [](const gph::HADRECAppModule& self, const int& bus_idx) -> py::object {
           int n(self.numLoads(bus_idx));
           return py::cast(n);
         })
#endif

#ifdef USE_OVERLOAD_CAST    
    .def("numStorage",
         py::overload_cast<>(&gph::HADRECAppModule::numStorage, py::const_))
    .def("numStorage",
         py::overload_cast<const int&>(&gph::HADRECAppModule::numStorage, py::const_))
#else
    .def("numStorage",
         [](const gph::HADRECAppModule& self) -> py::object {
           int n(self.numStorage());
           return py::cast(n);
         })
    .def("numStorage",
         [](const gph::HADRECAppModule& self, const int& bus_idx) -> py::object {
           int n(self.numStorage(bus_idx));
           return py::cast(n);
         })
#endif

#ifdef USE_OVERLOAD_CAST
    .def("numLines",
         py::overload_cast<>(&gph::HADRECAppModule::numLines, py::const_))
    .def("numLines",
         py::overload_cast<const int&>(&gph::HADRECAppModule::numLines, py::const_))
#else
    .def("numLines",
         [](const gph::HADRECAppModule& self) -> py::object {
           int n(self.numLines());
           return py::cast(n);
         })
    .def("numLines",
         [](const gph::HADRECAppModule& self, const int& branch_idx) -> py::object {
           int n(self.numLines(branch_idx));
           return py::cast(n);
         })

#endif
    
    .def("getBusInfoInt", 
         [](gph::HADRECAppModule& self, const int& bus_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           int result;
           bool ok = self.getBusInfo(bus_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBusInfoBool",
         [](gph::HADRECAppModule& self, const int& bus_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           bool result;
           bool ok = self.getBusInfo(bus_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBusInfoReal",
         [](gph::HADRECAppModule& self, const int& bus_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           double result;
           bool ok = self.getBusInfo(bus_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBusInfoString", 
         [](gph::HADRECAppModule& self, const int& bus_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           std::string result;
           bool ok = self.getBusInfo(bus_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBranchInfoInt", 
         [](gph::HADRECAppModule& self, const int& branch_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           int result;
           bool ok = self.getBranchInfo(branch_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBranchInfoBool",
         [](gph::HADRECAppModule& self, const int& branch_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           bool result;
           bool ok = self.getBranchInfo(branch_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBranchInfoReal",
         [](gph::HADRECAppModule& self, const int& branch_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           double result;
           bool ok = self.getBranchInfo(branch_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    .def("getBranchInfoString", 
         [](gph::HADRECAppModule& self, const int& branch_idx,
            const std::string& name, const std::optional<int> dev_idx) -> py::object {
           std::string result;
           bool ok = self.getBranchInfo(branch_idx, name, result, dev_idx.value());
           if (ok) {
             return py::cast(result);
           } else {
             return py::cast<py::none>(Py_None);
           }
         }, py::arg("bus_idx"), py::arg("name"), py::arg("dev_idx") = - 1)
    ;

  // These methods need to be reworked char * and/or optional args
  hadapp
    .def("scatterInjectionLoad",
         [](gph::HADRECAppModule& self, const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ) {
           self.scatterInjectionLoad(vbusNum, vloadP, vloadQ);
         },
         py::arg("vbusNum") = std::vector< int >(), py::arg("vloadP") = std::vector< double >(), py::arg("vloadQ") = std::vector< double >()
         )
    .def("scatterInjectionLoadNew",
         [](gph::HADRECAppModule& self, const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ) {
           self.scatterInjectionLoadNew(vbusNum, vloadP, vloadQ);
         },
         py::arg("vbusNum") = std::vector< int >(), py::arg("vloadP") = std::vector< double >(), py::arg("vloadQ") = std::vector< double >()
         )
    .def("scatterInjectionLoadNew_compensateY",
         [](gph::HADRECAppModule& self, const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ) {
           self.scatterInjectionLoadNew_compensateY(vbusNum, vloadP, vloadQ);
         },
         py::arg("vbusNum") = std::vector< int >(), py::arg("vloadP") = std::vector< double >(), py::arg("vloadQ") = std::vector< double >()
         )
    .def("scatterInjectionLoadNew_Norton",
         [](gph::HADRECAppModule& self, const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ, 
            const std::vector<double>& vimpedanceR, const std::vector<double>& vimpedanceI) {
           self.scatterInjectionLoadNew_Norton(vbusNum, vloadP, vloadQ, vimpedanceR, vimpedanceI);
         },
         py::arg("vbusNum") = std::vector< int >(), py::arg("vloadP") = std::vector< double >(), py::arg("vloadQ") = std::vector< double >(),
         py::arg("vimpedanceR") = std::vector< double >(), py::arg("vimpedanceI") = std::vector< double >()
         )
    .def("scatterInjectionLoadNewConstCur",
         [](gph::HADRECAppModule& self, const std::vector<int>& vbusNum, const std::vector<double>& vCurR, const std::vector<double>& vCurI) {
           self.scatterInjectionLoadNewConstCur(vbusNum, vCurR, vCurI);
         },
         py::arg("vbusNum") = std::vector< int >(), py::arg("vCurR") = std::vector< double >(), py::arg("vCurI") = std::vector< double >()
         )
    .def("initializeDynSimu",
         [](gph::HADRECAppModule& self, std::vector< gpds::Event > faults, int dscase_idx) {
           self.initializeDynSimu(faults, dscase_idx);
         },
         py::arg("faults") = std::vector< gpds::Event >(), py::arg("dscase_idx") = -1
         )
    .def("solvePowerFlowBeforeDynSimu",
         [](gph::HADRECAppModule& self, const std::string& s, int pfcase_idx) {
           self.solvePowerFlowBeforeDynSimu(s.c_str(), pfcase_idx);
         },
         py::arg("s") = "", py::arg("pfcase_idx") = -1
         )
    .def("exportPSSE23",
         [](gph::HADRECAppModule& self, const std::string& s) {
           self.exportPSSE23(s.c_str());
         },
         py::arg("s") = ""
         )
    .def("exportPSSE33",
         [](gph::HADRECAppModule& self, const std::string& s) {
           self.exportPSSE33(s.c_str());
         },
         py::arg("s") = ""
         )
    .def("exportPSSE34",
         [](gph::HADRECAppModule& self, const std::string& s) {
           self.exportPSSE34(s.c_str());
         },
         py::arg("s") = ""
         )
    .def("solvePowerFlowBeforeDynSimu_withFlag",
         [](gph::HADRECAppModule& self, const std::string& s, int pfcase_idx) {
           bool flag;
           flag = self.solvePowerFlowBeforeDynSimu_withFlag(s.c_str(), pfcase_idx);
           return flag;
         },
         py::arg("s") = "", py::arg("pfcase_idx") = -1
         )
    .def("readPowerFlowData",
         [](gph::HADRECAppModule& self, const std::string& s, int pfcase_idx) {
           self.readPowerFlowData(s.c_str(), pfcase_idx);
         },
         py::arg("s") = "", py::arg("pfcase_idx") = -1
         )
    .def("solvePowerFlow",
         [](gph::HADRECAppModule& self) {
           bool flag;
           flag = self.solvePowerFlow();
           return flag;
         }
         )
    .def("setWideAreaControlSignal",
         [](gph::HADRECAppModule& self, int bus_number, const std::string& genid, double wideAreaControlSignal) {
           self.setWideAreaControlSignal(bus_number, genid.c_str(), wideAreaControlSignal);
         },
         py::arg("bus_number") = "-1", py::arg("genid") = "", py::arg("wideAreaControlSignal") = 0.0
         )
    .def("modifyDataCollectionGenParam",
         [](gph::HADRECAppModule& self, int bus_no, const std::string& s, const std::string& par, double modvalue) {
           bool flag;
           flag = self.modifyDataCollectionGenParam(bus_no, s.c_str(), par.c_str(), modvalue);
           return flag;
         }
         )
    .def("modifyDataCollectionGenParam",
         [](gph::HADRECAppModule& self, int bus_no, const std::string& s, const std::string& par, int modvalue) {
           bool flag;
           flag = self.modifyDataCollectionGenParam(bus_no, s.c_str(), par.c_str(), modvalue);
           return flag;
         }
         )
    .def("modifyDataCollectionLoadParam",
         [](gph::HADRECAppModule& self, int bus_no, const std::string& s, const std::string& par, double modvalue) {
           bool flag;
           flag = self.modifyDataCollectionLoadParam(bus_no, s.c_str(), par.c_str(), modvalue);
           return flag;
         }
         )
    .def("modifyDataCollectionLoadParam",
         [](gph::HADRECAppModule& self, int bus_no, const std::string& s, const std::string& par, int modvalue) {
           bool flag;
           flag = self.modifyDataCollectionLoadParam(bus_no, s.c_str(), par.c_str(), modvalue);
           return flag;
         }
         )
    .def("modifyDataCollectionBusParam",
         [](gph::HADRECAppModule& self, int bus_no, const std::string& par, double modvalue) {
           bool flag;
           flag = self.modifyDataCollectionBusParam(bus_no, par.c_str(), modvalue);
           return flag;
         }
         )
    .def("modifyDataCollectionBusParam",
         [](gph::HADRECAppModule& self, int bus_no, const std::string& par, int modvalue) {
           bool flag;
           flag = self.modifyDataCollectionBusParam(bus_no, par.c_str(), modvalue);
           return flag;
         }
         )
    .def("modifyDataCollectionBranchParam",
         [](gph::HADRECAppModule& self, int bus1, int bus2, const std::string& s, const std::string& par, double modvalue) {
           bool flag;
           flag = self.modifyDataCollectionBranchParam(bus1, bus2, s.c_str(), par.c_str(), modvalue);
           return flag;
         }
         )		 
    .def("modifyDataCollectionBranchParam",
         [](gph::HADRECAppModule& self, int bus1, int bus2, const std::string& s, const std::string& par, int modvalue) {
           bool flag;
           flag = self.modifyDataCollectionBranchParam(bus1, bus2, s.c_str(), par.c_str(), modvalue);
           return flag;
         }
         )
    .def("fullInitializationBeforeDynSimuSteps",
         [](gph::HADRECAppModule& self, const std::string& s,
            const std::vector<gpds::Event>& BusFaults, int pfcase_idx, int dscase_idx) {
           self.fullInitializationBeforeDynSimuSteps(s.c_str(), BusFaults,
                                                     pfcase_idx, dscase_idx);
         },
         py::arg("s") = "",
         py::arg("BusFaults") = std::vector<gpds::Event>(),
         py::arg("pfcase_idx") = -1,
         py::arg("dscase_idx")
         )
    .def("setState",
         [](gph::HADRECAppModule& self, const int& bus_id, const std::string& dev_id,
            const std::string& device, const std::string& name,
            const double& value) {
           bool status;
           status = self.setState(bus_id, dev_id, device, name, value);
           return status;
         }
         )
    .def("getState",
         [](gph::HADRECAppModule& self, const int& bus_id, const std::string& dev_id,
            const std::string& device, const std::string& name) -> py::object {
           bool status;
           double value;
           status = self.getState(bus_id, dev_id, device, name, &value);

           // Value is returned if successful.
           if (status) {
             return  py::cast(value);
           }

           // Otherwise, None is returned. 
           return py::object();
         }
         )
    ;

  // This method returns a tuple containing 3 lists (int, string, int)

  hadapp
    .def("getObservationLists",
         [](gph::HADRECAppModule& self) {
           std::vector<int> obs_genBus;
           std::vector<std::string> obs_genIDs;
           std::vector<int> obs_loadBuses;
           std::vector<std::string> obs_loadIDs;
           std::vector<int> obs_busIDs;
           self.getObservationLists(obs_genBus, obs_genIDs,
                                    obs_loadBuses, obs_loadIDs, obs_busIDs);
           return py::make_tuple(obs_genBus, obs_genIDs,
                                 obs_loadBuses, obs_loadIDs, obs_busIDs);
         })
    .def("getObservationLists_withBusFreq",
         [](gph::HADRECAppModule& self) {
           std::vector<int> obs_genBus;
           std::vector<std::string> obs_genIDs;
           std::vector<int> obs_loadBuses;
           std::vector<std::string> obs_loadIDs;
           std::vector<int> obs_busIDs;
           std::vector<int> busfreqIDs;
           self.getObservationLists_withBusFreq(obs_genBus, obs_genIDs,
                                                obs_loadBuses, obs_loadIDs, obs_busIDs, busfreqIDs);
           return py::make_tuple(obs_genBus, obs_genIDs,
                                 obs_loadBuses, obs_loadIDs, obs_busIDs, busfreqIDs);
         })
    ;
             
  // These methods return a tuple on success or False on failure
  hadapp
    .def("getBusTotalLoadPower",
         [](gph::HADRECAppModule& self, const int& busid) -> py::object {
           double pg, qg;
           bool flag;
           flag = self.getBusTotalLoadPower(busid, pg, qg);
           if (flag) {
             return py::make_tuple(pg, qg);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getGeneratorPower",
         [](gph::HADRECAppModule& self, const int& busid, const std::string& genid) -> py::object {
           double pg, qg;
           bool flag;
           flag = self.getGeneratorPower(busid, genid, pg, qg);
           if (flag) {
             return py::make_tuple(pg, qg);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getPFSolutionSingleBus",
         [](gph::HADRECAppModule& self, const int& busid) -> py::object {
           double vmag, vangle;
           bool flag;
           flag = self.getPFSolutionSingleBus(busid, vmag, vangle);
           if (flag) {
             return py::make_tuple(vmag, vangle);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getDataCollectionGenParam",
         [](gph::HADRECAppModule& self, const int& busid, const std::string& genid, const std::string& genparam) -> py::object {
           double value;
           bool flag;
           flag = self.getDataCollectionGenParam(busid, genid, genparam, value);
           if (flag) {
             return py::float_(value);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getDataCollectionGenParam",
         [](gph::HADRECAppModule& self, const int& busid, const std::string& genid, const std::string& genparam) -> py::object {
           int value;
           bool flag;
           flag = self.getDataCollectionGenParam(busid, genid, genparam, value);
           if (flag) {
             return py::int_(value);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getDataCollectionLoadParam",
         [](gph::HADRECAppModule& self, const int& busid, const std::string& loadid, const std::string& loadparam) -> py::object {
           double value;
           bool flag;
           flag = self.getDataCollectionLoadParam(busid, loadid, loadparam, value);
           if (flag) {
             return py::float_(value);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getDataCollectionLoadParam",
         [](gph::HADRECAppModule& self, const int& busid, const std::string& loadid, const std::string& loadparam) -> py::object {
           int value;
           bool flag;
           flag = self.getDataCollectionLoadParam(busid, loadid, loadparam, value);
           if (flag) {
             return py::int_(value);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getDataCollectionBusParam",
         [](gph::HADRECAppModule& self, const int& busid, const std::string& param) -> py::object {
           double value;
           bool flag;
           flag = self.getDataCollectionBusParam(busid, param, value);
           if (flag) {
             return py::float_(value);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getDataCollectionBusParam",
         [](gph::HADRECAppModule& self, const int& busid, const std::string& param) -> py::object {
           int value;
           bool flag;
           flag = self.getDataCollectionBusParam(busid, param, value);
           if (flag) {
             return py::int_(value);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getDataCollectionBranchParam",
         [](gph::HADRECAppModule& self, const int& busid1, const int& busid2, const std::string& ckt, const std::string& param) -> py::object {
           double value;
           bool flag;
           flag = self.getDataCollectionBranchParam(busid1, busid2, ckt, param, value);
           if (flag) {
             return py::float_(value);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getDataCollectionBranchParam",
         [](gph::HADRECAppModule& self, const int& busid1, const int& busid2, const std::string& ckt, const std::string& param) -> py::object {
           int value;
           bool flag;
           flag = self.getDataCollectionBranchParam(busid1, busid2, ckt, param, value);
           if (flag) {
             return py::int_(value);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getZoneLoads",
         [](gph::HADRECAppModule& self) -> py::object {
           std::vector<double> load_p, load_q;
           std::vector<int> zone_id;
           bool flag = self.getZoneLoads(load_p, load_q, zone_id);
           if (flag) {
             return py::make_tuple(load_p, load_q, zone_id);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })
    .def("getZoneGeneratorPower",
         [](gph::HADRECAppModule& self) -> py::object {
           std::vector<double> generator_p, generator_q;
           std::vector<int> zone_id;
           bool flag = self.getZoneGeneratorPower(generator_p, generator_q, zone_id);
           if (flag) {
             return py::make_tuple(generator_p, generator_q, zone_id);
           } else {
             return py::cast<py::none>(Py_None);
           }
         })

    
    ;

}

