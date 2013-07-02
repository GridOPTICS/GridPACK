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
    int branch_idx, bus1_idx, bus2_idx;
    p_network->getBranchEndpoints(i, &idx1, &idx2);
    p_network->getBranch(i)->setBus1(p_network->getBus(idx1));
    p_network->getBranch(i)->setBus2(p_network->getBus(idx2));
    branch_idx = p_network->getGlobalBranchIndex(i); 
    bus1_idx = p_network->getGlobalBusIndex(idx1); 
    bus2_idx = p_network->getGlobalBusIndex(idx2); 
    p_network->getBranch(i)->setGlobalIndices(branch_idx, bus1_idx, bus2_idx); 
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
    int bus_idx;
    bus_idx = p_network->getGlobalBusIndex(i); 
    p_network->getBus(i)->setGlobalIndex(bus_idx);
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

/**
 * Set up the exchange buffers so that they work correctly. This should only
 * be called after the network topology has been specified
 */
void gridpack::factory::BaseFactory::setExchange(void)
{
  int busXCSize, branchXCSize;
  int nbus, nbranch;

  nbus = p_network->numBuses();
  nbranch = p_network->numBranches();

  // Get size of bus and branch exchange buffers from first local bus and branch
  // components. These must be the same for all bus and branch components

  if (nbus > 0) {
    busXCSize = p_network->getBus(0)->getXCBufSize();
  } else {
    busXCSize = 0;
  }
  p_network->allocXCBus(busXCSize);

  if (nbranch > 0){
    branchXCSize = p_network->getBranch(0)->getXCBufSize();
  }
  p_network->allocXCBranch(branchXCSize);

  // Buffers have been allocated in network. Now associate buffers from network
  // back to individual components
  int i;
  for (i=0; i<nbus; i++) {
    p_network->getBus(i)->setXCBuf(p_network->getXCBusBuffer(i));
  }
  for (i=0; i<nbranch; i++) {
    p_network->getBranch(i)->setXCBuf(p_network->getXCBranchBuffer(i));
  }
}

/**
 * Set the mode for all BaseComponent objects in the network.
 * @param mode: integer representing desired mode
 */
void gridpack::factory::BaseFactory::setMode(int mode)
{
  int nbus, nbranch;

  nbus = p_network->numBuses();
  nbranch = p_network->numBranches();

  int i;
  for (i=0; i<nbus; i++) {
    p_network->getBus(i)->setMode(mode);
  }
  for (i=0; i<nbranch; i++) {
    p_network->getBranch(i)->setMode(mode);
  }
}
