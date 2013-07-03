// -------------------------------------------------------------
/**
 * @file   pf_factory.cpp
 * @author Bruce Palmer
 * @date   July 1, 2013
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
#include "gridpack/applications/dynsim/pf_components.hpp"
#include "gridpack/applications/dynsim/pf_factory.hpp"

// Powerflow factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
gridpack::dynsim::DynSimFactory::DynSimFactory(
         boost::shared_ptr<gridpack::network::BaseNetwork
         <gridpack::component::BaseBusComponent,
         gridpack::component::BaseBranchComponent> > network)
         : gridpack::factory::BaseFactory(network)
{
  p_network = network;
}

/**
 * Basic destructor
 */
gridpack::dynsim::DynSimFactory::~DynSimFactory()
{
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::dynsim::DynSimFactory::setYBus(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // Invoke setYBus method on all bus objects
  for (i=0; i<numBus; i++) {
    (dynamic_cast<gridpack::dynsim::DynSimBus*>(p_network->getBus(i).get()))->setYBus();
  }

  // Invoke setYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
    (dynamic_cast<gridpack::dynsim::DynSimBranch*>(p_network->getBranch(i).get()))->setYBus();
  }
}
