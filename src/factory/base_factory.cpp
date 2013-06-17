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
#include "gridpack/factory/base_factory.hpp"

// Base factory class that contains functions that are generic to all
// applications.

/**
 * Constructor
 */
gridpack::factory::BaseFactory::BaseFactory(
  boost::shared_ptr<gridpack::network::BaseNetwork<gridpack::component::BaseBusComponent,
                    gridpack::component::BaseBranchComponent> > network)
{
  p_network = network;
}

/**
 * Destructor
 */
gridpack::factory::BaseFactory::~BaseFactory(void)
{
}

/**
 * Set pointers in each bus and branch component so that it points to
 * connected buses and branches. This routine operates on the generic
 * BaseBusComponent and BaseBranchComponent interfaces.
 */
void gridpack::factory::BaseFactory::setComponents(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i, j;

  // Set pointers for buses at either end of each branch
  for (i=0; i<numBranch; i++) {
    int idx1, idx2;
    p_network->getBranchEndpoints(i, &idx1, &idx2);
    p_network->getBranch(i)->setBus1(p_network->getBus(idx1));
    p_network->getBranch(i)->setBus2(p_network->getBus(idx2));
  }

  // Set pointers for branches and buses connected to each bus
  for (i=0; i<numBus; i++) {
    p_network->getBus(i)->clearBuses();
    std::vector<int> nghbrBus = p_network->getConnectedBuses(i);
    for (j=0; j<nghbrBus.size(); j++) {
      p_network->getBus(i)->addBus(p_network->getBus(nghbrBus[j]));
    }
    p_network->getBus(i)->clearBranches();
    std::vector<int> nghbrBranch = p_network->getConnectedBranches(i);
    for (j=0; j<nghbrBranch.size(); j++) {
      p_network->getBus(i)->addBranch(p_network->getBranch(nghbrBranch[j]));
    }
  }
}

/**
 * Generic method that invokes the "load" method on all branches and buses
 * to move data from the DataCollection objects on the network into the
 * corresponding buses and branches
 */
void gridpack::factory::BaseFactory::load(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // Invoke load method on all bus objects
  for (i=0; i<numBus; i++) {
    p_network->getBus(i)->load(p_network->getBusData(i));
  }

  // Invoke load method on all branch objects
  for (i=0; i<numBranch; i++) {
    p_network->getBranch(i)->load(p_network->getBranchData(i));
  }
}
