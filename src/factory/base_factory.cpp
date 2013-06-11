// -------------------------------------------------------------
/**
 * @file   base_factory.cpp
 * @author Bruce Palmer
 * @date   June 3, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"

// Base factory class that contains functions that are generic to all
// applications.

/**
 * Constructor
 */
BaseFactory(void)
{
}

/**
 * Destructor
 */
virtual ~BaseFactory(void)
{
}

/**
 * Set pointers in each bus and branch component so that it points to
 * connected buses and branches. This routine operates on the generic
 * BaseBusComponent and BaseBranchComponent interfaces.
 * @param network: The network that contains the components that need to be
 * set
 */
virtual void setComponents(gridpack::network::BaseNetwork
                               <gridpack::component::BaseBusComponent>,
                               gridpack::component::BaseBranchComponent> >
                               *network)
{
  int numBus = network->numBuses();
  int numBranch = network->numBranches();
  int i, j;

  // Set pointers for buses at either end of each branch
  for (i=0; i<numBranch; i++) {
    int idx1, idx2;
    network->getBranchEndpoints(i, &idx1, &idx2);
    network->getBranch(i)->setBus1(network->getBus(idx1));
    network->getBranch(i)->setBus2(network->getBus(idx2));
  }

  // Set pointers for branches and buses connected to each bus
  for (i=0; i<numBus; i++) {
    network->getBus(i)->clearBuses();
    std::vector<int> nghbrBus = network->getConnectedBuses(i);
    for (j=0; j<nghbrBus.size(); j++) {
      network->getBus(i)->addBus(network->getBus(nghbrBus[j]));
    }
    network->getBus(i)->clearBranches();
    std::vector<int> nghbrBranch = network->getConnectedBranches(i);
    for (j=0; j<nghbrBranch.size(); j++) {
      network->getBus(i)->addBranch(network->getBranch(nghbrBranch[j]));
    }
  }
}
