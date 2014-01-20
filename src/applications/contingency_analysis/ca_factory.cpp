/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ca_factory.cpp
 * @author Yousu Chen 
 * @date   January 20, 2014
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
#include "gridpack/applications/contingency_analysis/ca_components.hpp"
#include "gridpack/applications/contingency_analysis/ca_factory.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"

namespace gridpack {
namespace contingency_analysis {

// Contingency analysis factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
CAFactory::CAFactory(CAFactory::NetworkPtr network)
  : gridpack::factory::BaseFactory<CANetwork>(network)
{
  p_network = network;
}

/**
 * Basic destructor
 */
gridpack::contingency_analysis::CAFactory::~CAFactory()
{
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::contingency_analysis::CAFactory::setYBus(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // Invoke setYBus method on all bus objects
  for (i=0; i<numBus; i++) {
    (dynamic_cast<CABus*>(p_network->getBus(i).get()))->setYBus();
  }

  // Invoke setYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
    (dynamic_cast<CABranch*>(p_network->getBranch(i).get()))->setYBus();
  }
}

} // namespace contingency_analysis 
} // namespace gridpack
