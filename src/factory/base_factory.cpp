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
 * BaseBusComponent and BaseBranchComponent interfaces.It also sets some
 * indices in MatVecInterface for each component.
 */
void gridpack::factory::BaseFactory::setComponents(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i, j;
  int idx1, idx2;

  // Set pointers for buses at either end of each branch
  int numActiveBranch = 0;
  for (i=0; i<numBranch; i++) {
    int branch_idx, bus1_idx, bus2_idx;
    p_network->getBranchEndpoints(i, &idx1, &idx2);
    p_network->getBranch(i)->setBus1(p_network->getBus(idx1));
    p_network->getBranch(i)->setBus2(p_network->getBus(idx2));
    branch_idx = p_network->getGlobalBranchIndex(i); 
    bus1_idx = p_network->getGlobalBusIndex(idx1); 
    bus2_idx = p_network->getGlobalBusIndex(idx2); 
    p_network->getBranch(i)->setMatVecIndices(bus1_idx, bus2_idx); 
    if (p_network->getActiveBranch(i)) numActiveBranch++;
  }

  // Set pointers for branches and buses connected to each bus
  int numActiveBus = 0;
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
    p_network->getBus(i)->setMatVecIndex(bus_idx);
    if (p_network->getActiveBranch(i)) numActiveBus++;
  }

  // Come up with a set of global indices for each component so that the buses
  // and branches are consecutively numbered on each component and the indices
  // for all active components are unique. These are used in the mapper routine.
  // First get the number of active buses and branches on each process and
  // broadcast this to all other processes.
  int nprocs = GA_Nnodes();
  int me = GA_Nodeid();
  int *activeBus = new int[nprocs];
  int *activeBranch = new int[nprocs];

  for (i=0; i<nprocs; i++) {
    activeBus[i] = 0;
    activeBranch[i] = 0;
  }

  activeBus[me] = numActiveBus;
  activeBranch[me] = numActiveBranch;
  GA_Igop(activeBus,nprocs,"+");
  GA_Igop(activeBranch,nprocs,"+");

  // Create indices for buses. Start by creating a global array with an entry
  // for each bus
  int one = 1;
  int ntot = 0;
  int offset = 0;
  for (i=0; i<nprocs; i++) {
    ntot += activeBus[i];
    if (i<me) offset += activeBus[i];
  }
  int g_bus = GA_Create_handle();
  GA_Set_data(g_bus, one, &ntot, C_INT);
  if (!GA_Allocate(g_bus)) {
    // TODO: some kind of error
  }
  GA_Zero(g_bus);

  // Scatter new index values into the global index locations for the buses
  int *ibus_val = new int[numActiveBus];
  int **ibus_idx = new int*[numActiveBus];
  int icnt = 0;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      p_network->getBus(i)->setMatVecIndex(icnt+offset);
      ibus_idx[icnt] = new int;
      *(ibus_idx[icnt]) = p_network->getGlobalBusIndex(i);
      ibus_val[icnt] = offset+icnt;
      icnt++;
    }
  }
  NGA_Scatter(g_bus,ibus_val,ibus_idx,numActiveBus);
  GA_Sync();
  for (i=0; i<numActiveBus; i++) {
    delete ibus_idx[i];
  }
  delete [] ibus_idx;
  delete [] ibus_val;

  // Now gather values for both active and inactive buses
  ibus_val = new int[numBus];
  ibus_idx = new int*[numBus];
  for (i=0; i<numBus; i++) {
    ibus_idx[i] = new int;
    *(ibus_idx[icnt]) = p_network->getGlobalBusIndex(i);
  }
  NGA_Gather(g_bus,ibus_val,ibus_idx,numBus);
  // Assign the MatVecIndex for the bus and clean up arrays
  for (i=0; i<numBus; i++) {
    p_network->getBus(i)->setMatVecIndex(ibus_val[i]);
    delete ibus_idx[i];
  }
  delete [] ibus_idx;
  delete [] ibus_val;
  delete [] activeBus;
  delete [] activeBranch;
  GA_Destroy(g_bus);

  // Finish by assigning MatVecIndices for the branches
  for (i=0; i<numBranch; i++) {
    p_network->getBranch(i)->getBus1()->getMatVecIndex(&idx1);
    p_network->getBranch(i)->getBus2()->getMatVecIndex(&idx2);
    p_network->getBranch(i)->setMatVecIndices(idx1,idx2);
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
