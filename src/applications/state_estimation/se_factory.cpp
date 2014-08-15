/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_factory.cpp
 * @author Yousu Chen 
 * @date   2/24/2014 
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "se_factory.hpp"

namespace gridpack {
namespace state_estimation {

// State estimation factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
SEFactory::SEFactory(SEFactory::NetworkPtr network)
  : gridpack::factory::BaseFactory<SENetwork>(network)
{
  p_network = network;
}

/**
 * Basic destructor
 */
gridpack::state_estimation::SEFactory::~SEFactory()
{
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::state_estimation::SEFactory::setYBus(void)
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
void gridpack::state_estimation::SEFactory::setMeasurements(
    std::vector<gridpack::state_estimation::Measurement> measurements)
{
  int me = p_network->communicator().rank();
  int size = measurements.size();
  std::vector<int> bus_keys;
  std::vector<int> branch_ids;
  std::vector<std::pair<int,int> > branch_keys;
  std::vector<gridpack::state_estimation::Measurement> bus_meas;
  std::vector<gridpack::state_estimation::Measurement> branch_meas;
  int i, nbus, nbranch, idx1, idx2;
  std::pair<int,int> key;
  for (i=0; i<size; i++) {
    std::string meas_type = measurements[i].p_type;
    if (meas_type == "VM" || meas_type == "PI" ||
        meas_type == "PJ" || meas_type == "QI" ||
        meas_type == "QJ" || meas_type == "VA") {
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
}

} // namespace state_estimation
} // namespace gridpack
