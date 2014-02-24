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
#include "gridpack/network/base_network.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/applications/state_estimation/se_components.hpp"
#include "gridpack/applications/state_estimation/se_factory.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"

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

} // namespace state_estimation
} // namespace gridpack
