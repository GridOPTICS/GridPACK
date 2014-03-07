/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_app.cpp
 * @author Shuangshuang Jin
 * @date   March 06, 2014
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#include "gridpack/math/matrix.hpp"
#include "gridpack/math/vector.hpp"
#include "gridpack/math/linear_solver.hpp"
#include "gridpack/math/linear_matrix_solver.hpp"
#include "gridpack/applications/dynamic_simulation_2/ds2_app.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/applications/dynamic_simulation_2/ds2_factory.hpp"

// Calling program for dynamic simulation application

/**
 * Basic constructor
 */
gridpack::dynamic_simulation::DSApp::DSApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::DSApp::~DSApp(void)
{
}

/**
 * Execute application
 */
void gridpack::dynamic_simulation::DSApp::execute(int argc, char** argv)
{
  gridpack::parallel::Communicator world;
  boost::shared_ptr<DSNetwork> network(new DSNetwork(world));

  // read configuration file
  gridpack::utility::Configuration *config = gridpack::utility::Configuration::configuration();
  if (argc >= 2 && argv[1] != NULL) {
    char inputfile[256];
    sprintf(inputfile,"%s",argv[1]);
    config->open(inputfile,world);
  } else {
    config->open("input.xml",world);
  }
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Dynamic_simulation");
  std::string filename = cursor->get("networkConfiguration",
      "No network configuration specified");
  double sim_time = cursor->get("simulationTime",0.0);
  if (sim_time == 0.0) {
    // TODO: some kind of error
  }
  double time_step = cursor->get("timeStep",0.0);
  if (time_step == 0.0) {
    // TODO: some kind of error
  }

  // Read in information about fault events and store them in internal data
  // structure
  cursor = config->getCursor("Configuration.Dynamic_simulation.faultEvents");
  gridpack::utility::Configuration::ChildCursors events;
  if (cursor) cursor->children(events);
  std::vector<gridpack::dynamic_simulation::DSBranch::Event>
     faults = setFaultEvents(events,network); 

  // load input file
  gridpack::parser::PTI23_parser<DSNetwork> parser(network);
  parser.parse(filename.c_str());

  // partition network
  network->partition();

  // Create serial IO object to export data from buses or branches
  gridpack::serial_io::SerialBusIO<DSNetwork> busIO(256, network);
  gridpack::serial_io::SerialBranchIO<DSNetwork> branchIO(128, network);
  char ioBuf[128];

  // create factory
  gridpack::dynamic_simulation::DSFactory factory(network);
  factory.load();

  // set network components using factory
  factory.setComponents();

  // set YBus components so that you can create Y matrix  
  factory.setYBus();

  if (!factory.checkGen()) {
    busIO.header("Missing generators on at least one processor\n");
    return;
  }

  factory.setMode(YBUS);
  gridpack::mapper::FullMatrixMap<DSNetwork> ybusMap(network);
  boost::shared_ptr<gridpack::math::Matrix> orgYbus = ybusMap.mapToMatrix();
  branchIO.header("\n=== orginal ybus: ============\n");
  orgYbus->print();

  // Form constant impedance load admittance yl for all buses and add it to 
  // system Y matrix: ybus = ybus + yl
  factory.setMode(YL);
  boost::shared_ptr<gridpack::math::Matrix> ybus = ybusMap.mapToMatrix();
  branchIO.header("\n=== ybus after added yl: ============\n");
  ybus->print();
}

/**
 * Utility function to convert faults that are in event list into
 * internal data structure that can be used by code
 * @param cursors list of cursors pointing to individual events in input
 * deck
 * @return list of event data structures
 */
std::vector<gridpack::dynamic_simulation::DSBranch::Event>
   gridpack::dynamic_simulation::DSApp::setFaultEvents(
   std::vector<gridpack::utility::Configuration::CursorPtr > events,
   boost::shared_ptr<DSNetwork> network)
{
  int size = events.size();
  int idx;
  std::vector<gridpack::dynamic_simulation::DSBranch::Event> faults;
  // Parse fault events
  for (idx=0; idx<size; idx++) {
    gridpack::dynamic_simulation::DSBranch::Event event;
    event.start = events[idx]->get("beginFault",0.0);
    event.end = events[idx]->get("endFault",0.0);
    std::string indices = events[idx]->get("faultBranch","0 0");
    //Parse indices to get from and to indices of branch
    int ntok1 = indices.find_first_not_of(' ',0);
    int ntok2 = indices.find(' ',ntok1);
    if (ntok2 - ntok1 > 0 && ntok1 != std::string::npos && ntok2 !=
        std::string::npos) {
      event.from_idx = atoi(indices.substr(ntok1,ntok2-ntok1).c_str());
      ntok1 = indices.find_first_not_of(' ',ntok2);
      ntok2 = indices.find(' ',ntok1);
      if (ntok1 != std::string::npos && ntok1 < indices.length()) {
        if (ntok2 == std::string::npos) {
          ntok2 = indices.length();
        }
        event.to_idx = atoi(indices.substr(ntok1,ntok2-ntok1).c_str());
      } else {
        event.from_idx = 0;
        event.to_idx = 0;
      }
    } else {
      event.from_idx = 0;
      event.to_idx = 0;
    }
    event.step = events[idx]->get("timeStep",0.0);
    if (event.step != 0.0 && event.end != 0.0 && event.from_idx != event.to_idx) {
      faults.push_back(event);
    }
  }
#if 0
  // Find local indices of branches on which faults occur. Start by constructing
  // a map object that maps to local branch indices using branch index pairs as the key
  std::map<std::pair<int, int>, int> pairMap;
  int numBranch = network->numBranches();
  for (idx = 0; idx<numBranch; idx++) {
    // Only set branch index for locally held branches
    if (network->getActiveBranch(idx)) {
      int idx1, idx2;
      network->getOriginalBranchEndpoints(idx, &idx1, &idx2);
      std::pair<int, int> branch_pair(idx1, idx2);
      pairMap.insert(std::pair<std::pair<int, int>, int>(branch_pair,idx));
    }
  }
  // run through all events and see if the branch exists on this processor. If
  // it does, then set the branch_idx member to the local branch index.
  size = faults.size();
  for (idx=0; idx<size; idx++) {
    std::map<std::pair<int, int>, int>::iterator it;
    std::pair<int, int> branch_pair(faults[idx].from_idx, faults[idx].to_idx);
    it = pairMap.find(branch_pair);
    if (it != pairMap.end()) {
      faults[idx].branch_idx = it->second;
    }
  }
#endif
  return faults;
}
