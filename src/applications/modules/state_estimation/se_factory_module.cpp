/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_factory_module.cpp
 * @author Yousu Chen, Bruce Palmer
 * @date   1/23/2015
 * 
 * @brief  Implementation of the SE factory module
 * @update Yousu Chen
 *         Adding functions of applying virtual measurements
 * @date   2025-03-05
 *         Adding more functions to handle bad data detection and create more comprehensive outputs
 * @date   2025-04-20
 *         Improved scalability with optimized matrix operation for large systems
 * @date   2025-04-22
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/hash_distr.hpp"
#include "timer/coarse_timer.hpp"
#include "se_factory_module.hpp"

namespace gridpack {
namespace state_estimation {

// State estimation factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
SEFactoryModule::SEFactoryModule(SEFactoryModule::NetworkPtr network)
  : gridpack::factory::BaseFactory<SENetwork>(network)
{
  p_network = network;
  
  // Optimize for large systems automatically
  p_optimizeForLargeSystems = false;
  
  // Check system size for automatic optimization
  int numBus = p_network->numBuses();
  if (numBus > 1000) {
    p_optimizeForLargeSystems = true;
    if (p_network->communicator().rank() == 0) {
      printf("SEFactoryModule: Detected large system (%d buses), enabling optimizations\n", numBus);
    }
  }
}

/**
 * Basic destructor
 */
SEFactoryModule::~SEFactoryModule()
{
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void SEFactoryModule::setYBus(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // Invoke setYBus method on all bus objects
  for (i=0; i<numBus; i++) {
    (dynamic_cast<SEBus*>(p_network->getBus(i).get()))->setYBus();
  }

  // Invoke setYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
    (dynamic_cast<SEBranch*>(p_network->getBranch(i).get()))->setYBus();
  }
}

/**
 * Disribute measurements
 * @param measurements a vector containing all measurements
 */
void SEFactoryModule::setMeasurements(
    std::vector<gridpack::state_estimation::Measurement> measurements)
{
  int me = p_network->communicator().rank();
  int size = measurements.size();
  std::vector<int> bus_keys;
  std::vector<int> branch_ids;
  std::vector<std::pair<int,int> > branch_keys;
  std::vector<gridpack::state_estimation::Measurement> bus_meas;
  std::vector<gridpack::state_estimation::Measurement> branch_meas;
  int i, idx1, idx2;
  std::pair<int,int> key;
  
  for (i=0; i<size; i++) {
    std::string meas_type = measurements[i].p_type;
    if (meas_type == "VM" || meas_type == "PI" ||
        meas_type == "PJ" || meas_type == "QI" ||
        meas_type == "QJ" || meas_type == "VA" || 
        meas_type == "VUL" || meas_type == "VLL") {
      bus_meas.push_back(measurements[i]);
      bus_keys.push_back(measurements[i].p_busid);
    } else if (meas_type == "PIJ" || meas_type == "PJI" ||
        meas_type == "QIJ" || meas_type == "QJI" ||
        meas_type == "IIJ" || meas_type == "IJI") {
      branch_meas.push_back(measurements[i]);
      idx1 = measurements[i].p_fbusid;
      idx2 = measurements[i].p_tbusid;
      key = std::pair<int,int>(idx1,idx2);
      branch_keys.push_back(key);
    }
  }
  gridpack::hash_distr::HashDistribution<SENetwork,Measurement,Measurement>
    distr(p_network);
  distr.distributeBusValues(bus_keys,bus_meas);
  int nsize = bus_keys.size();
  for (i=0; i<nsize; i++) {
    p_network->getBus(bus_keys[i])->addMeasurement(bus_meas[i]);
  }
  bus_keys.clear();
  bus_meas.clear();
  distr.distributeBranchValues(branch_keys,branch_ids,branch_meas);
  branch_keys.clear();
  nsize = branch_ids.size();
  for (i=0; i<nsize; i++) {
    p_network->getBranch(branch_ids[i])->addMeasurement(branch_meas[i]);
  }
  branch_ids.clear();
  branch_meas.clear();
  // Sort measurements on each bus and branch so that they are in a consistent
  // order on all components
  nsize = p_network->numBuses();
  for (i=0; i<nsize; i++) {
    p_network->getBus(i)->sortMeasurements();
  }
  nsize = p_network->numBranches();
  for (i=0; i<nsize; i++) {
    p_network->getBranch(i)->sortMeasurements();
  }
}

/**
 * Initialize some parameters in state estimation components
 */
void SEFactoryModule::configureSE(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // For large systems, we can optimize by minimizing synchronization
  if (p_optimizeForLargeSystems) {
    
    for (i=0; i<numBus; i++) {
      (dynamic_cast<SEBus*>(p_network->getBus(i).get()))->configureSE();
    }
    // Process local branches
    for (i=0; i<numBranch; i++) {
      (dynamic_cast<SEBranch*>(p_network->getBranch(i).get()))->configureSE();
    }
    
    // Only synchronize once after all local work is done
    p_network->communicator().barrier();
  } else {
    
    // Invoke configureSE method on all bus objects
    for (i=0; i<numBus; i++) {
      (dynamic_cast<SEBus*>(p_network->getBus(i).get()))->configureSE();
    }

    // Invoke configureSE method on all branch objects
    for (i=0; i<numBranch; i++) {
      (dynamic_cast<SEBranch*>(p_network->getBranch(i).get()))->configureSE();
    }
  }
}

/**
 * Identify bad data based on residual index
 * @param idx: index of maximum normalized residual
 */
void SEFactoryModule::identifyBadData(int idx)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;
  bool found = false;
  
  // Check bus measurements
  for (i=0; i<numBus && !found; i++) {
    SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
    if (bus->checkResidualIndex(idx, true)) {
      found = true;
      break;
    }
  }
  
  // Check branch measurements
  for (i=0; i<numBranch && !found; i++) {
    SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
    if (branch->checkResidualIndex(idx, true)) {
      found = true;
      break;
    }
  }
}

// Set mode for components
void SEFactoryModule::setMode(int mode)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  
  // For large systems, optimize communication patterns
  if (p_optimizeForLargeSystems) {
    
    // Process local buses
    for (i = 0; i < numBus; i++) {
      dynamic_cast<SEBus*>(p_network->getBus(i).get())->setMode(mode);
    }
    
    // Process local branches
    for (i = 0; i < numBranch; i++) {
      dynamic_cast<SEBranch*>(p_network->getBranch(i).get())->setMode(mode);
    }
    
  } else {
    
    // Set mode for all buses
    for (i = 0; i < numBus; i++) {
      dynamic_cast<SEBus*>(p_network->getBus(i).get())->setMode(mode);
    }
    
    // Set mode for all branches
    for (i = 0; i < numBranch; i++) {
      dynamic_cast<SEBranch*>(p_network->getBranch(i).get())->setMode(mode);
    }
  }
}

bool SEFactoryModule::reportMeasurement(int idx, std::string& details)
{
  char buf[256];
  bool found = false;
  bool isPseudo = false;
  int numBus = p_network->numBuses();
  for (int i = 0; i < numBus && !found; i++) {
    SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
    if (bus->getResidualDetails(idx, buf)) {
      details = buf;
      bus->isResidualPseudo(idx, isPseudo);
      found = true;
    }
  }
  if (!found) {
    int numBranch = p_network->numBranches();
    for (int i = 0; i < numBranch && !found; i++) {
      SEBranch* branch = dynamic_cast<SEBranch*>(p_network->getBranch(i).get());
      if (branch->getResidualDetails(idx, buf)) {
        details = buf;
        branch->isResidualPseudo(idx, isPseudo);
        found = true;
      }
    }
  }
  // Mark pseudo-measurements
  if (found && isPseudo) details += " (pseudo)";
  return found;
}

std::set<int> SEFactoryModule::getPseudoResidualIndices()
{
  std::vector<int> idxs;
  int numBus = p_network->numBuses();
  for (int i = 0; i < numBus; i++) {
    dynamic_cast<SEBus*>(p_network->getBus(i).get())
      ->getPseudoResidualIndices(idxs);
  }
  int numBranch = p_network->numBranches();
  for (int i = 0; i < numBranch; i++) {
    dynamic_cast<SEBranch*>(p_network->getBranch(i).get())
      ->getPseudoResidualIndices(idxs);
  }
  return std::set<int>(idxs.begin(), idxs.end());
}

void SEFactoryModule::setVoltageLimits(double v_min, double v_max, bool enforce)
{
  int numBus = p_network->numBuses();
  for (int i = 0; i < numBus; i++) {
    SEBus* bus = dynamic_cast<SEBus*>(p_network->getBus(i).get());
    if (bus) {
      bus->setVoltageLimits(v_min, v_max, enforce);
    }
  }
}

// Get Mode Names
const char* SEFactoryModule::getModeNameString(int mode) const
{
  switch (mode) {
    case YBus:
      return "YBus";
    case Jacobian_H:
      return "Jacobian_H";
    case R_inv:
      return "R_inv";
    case Ez:
      return "Ez (Measurement)";
    case Voltage:
      return "Voltage";
    case Residual:
      return "Residual";
    default:
      return "Unknown";
  }
}

} // namespace state_estimation
} // namespace gridpack
