/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   network_c.cpp
 * @author Bruce Palmer
 * @date   2014-08-4 11:05:08 d3g293
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
#include "gridpack/parallel/communicator.hpp"
#include "gridpack/network/base_network.hpp"
#include "../component/fortran_component.hpp"
//#include "fortran_network.hpp"

typedef gridpack::network::BaseNetwork<
    gridpack::fortran_component::FortranBusComponent,
    gridpack::fortran_component::FortranBranchComponent>
    FortranNetwork;

struct networkWrapper {
  boost::shared_ptr<FortranNetwork> network;
};

typedef gridpack::fortran_component::FortranBusComponent FortranBus;

/**
 * Create a new network
 * @param network pointer to new Fortran network object
 * @param comm communicator
 */
extern "C" void network_create(networkWrapper **wnetwork,
    gridpack::parallel::Communicator *comm)
{
  *wnetwork = new networkWrapper;
  (*wnetwork)->network.reset(new FortranNetwork(*comm));
}

/**
 * Destroy old network
 * @param network pointer to Fortran network object
 */
extern "C" void network_destroy(networkWrapper **wnetwork)
{
  delete (*wnetwork);
  *wnetwork = NULL;
}

/**
 * Add a bus locally to the network
 * @param network pointer to Fortran network object
 * @param idx original index of bus
 */
extern "C" void network_add_bus(networkWrapper *wnetwork, int idx)
{
  wnetwork->network->addBus(idx);
}

/**
 * Add a branch locally to the network. A branch is defined by
 * buses at either end
 * @param network pointer to Fortran network object
 * @param idx1 original bus index of "from" bus
 * @param idx2 original bus index of "to" bus
 */
extern "C" void network_add_branch(networkWrapper *wnetwork, int idx1, int idx2)
{
  wnetwork->network->addBranch(idx1,idx2);
}

/**
 * Number of local buses (both active and inactive) on processor
 * @param network pointer to Fortran network object
 * @return number of buses
 */
extern "C" int network_num_buses(networkWrapper *wnetwork)
{
  return wnetwork->network->numBuses();
}

/**
 * Return the total number of buses in the entire network
 * @param network pointer to Fortran network object
 * @return total number of buses
 */
extern "C" int network_total_buses(networkWrapper *wnetwork)
{
  return wnetwork->network->totalBuses();
}

/**
 * Number of local branches (both active and inactive) on processor
 * @param network pointer to Fortran network object
 * @return number of branches
 */
extern "C" int network_num_branches(networkWrapper *wnetwork)
{
  return wnetwork->network->numBranches();
}

/**
 * Return the total number of branches in the entire network
 * @param network pointer to Fortran network object
 * @return total number of branches
 */
extern "C" int network_total_branches(networkWrapper *wnetwork)
{
  return wnetwork->network->totalBranches();
}

/**
 * Designate a bus as a reference bus.
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 */
extern "C" void network_set_reference_bus(networkWrapper *wnetwork, int idx)
{
  idx--;
  wnetwork->network->setReferenceBus(idx);
}

/**
 * Return index of reference bus.
 * @param network pointer to Fortran network object
 * @return local index of reference bus. If reference bus is not on this
 * processor then return -1.
 */
extern "C" int network_get_reference_bus(networkWrapper *wnetwork)
{
  int idx = wnetwork->network->getReferenceBus();
  if (idx >= 0) idx++;
  return idx;
}

/**
 * Set the original index of the bus (from configuration file)
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @param o_idx original index assigned to bus
 * @return false if no bus exists for idx
 */
extern "C" bool network_set_original_bus_index(networkWrapper *wnetwork, int idx, int o_idx)
{
  idx--;
  return wnetwork->network->setOriginalBusIndex(idx,o_idx);
}

/**
 * Set the global index of the bus
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @param g_idx global index to be assigned to bus
 * @return false if no bus exists for idx
 */
extern "C" bool network_set_global_bus_index(networkWrapper *wnetwork, int idx, int g_idx)
{
  idx--;
  g_idx--;
  return wnetwork->network->setGlobalBusIndex(idx,g_idx);
}

/**
 * Set the global index of the branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param g_idx global index to be assigned to branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_global_branch_index(networkWrapper *wnetwork, int idx, int g_idx)
{
  idx--;
  g_idx--;
  return wnetwork->network->setGlobalBranchIndex(idx,g_idx);
}

/**
 * Set the original index of the bus at the "from" end of branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param b_idx original index of "from" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_original_bus_index1(networkWrapper *wnetwork, int idx, int b_idx)
{
  idx--;
  return wnetwork->network->setOriginalBusIndex1(idx,b_idx);
}

/**
 * Set the original index of the bus at the "to" end of branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param b_idx original index of "to" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_original_bus_index2(networkWrapper *wnetwork, int idx, int b_idx)
{
  idx--;
  return wnetwork->network->setOriginalBusIndex2(idx,b_idx);
}

/**
 * Set the global index of the bus at the "from" end of branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param b_idx global index of "from" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_global_bus_index1(networkWrapper *wnetwork, int idx, int b_idx)
{
  idx--;
  b_idx--;
  return wnetwork->network->setGlobalBusIndex1(idx,b_idx);
}

/**
 * Set the global index of the bus at the "to" end of branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param b_idx global index of "to" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_global_bus_index2(networkWrapper *wnetwork, int idx, int b_idx)
{
  idx--;
  b_idx--;
  return wnetwork->network->setGlobalBusIndex2(idx,b_idx);
}

/**
 * Set the local index of the bus at the "from" end of branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param b_idx local index of "from" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_local_bus_index1(networkWrapper *wnetwork, int idx, int b_idx)
{
  idx--;
  b_idx--;
  return wnetwork->network->setLocalBusIndex1(idx,b_idx);
}

/**
 * Set the local index of the bus at the "to" end of branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param b_idx local index of "to" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_local_bus_index2(networkWrapper *wnetwork, int idx, int b_idx)
{
  idx--;
  b_idx--;
  return wnetwork->network->setLocalBusIndex2(idx,b_idx);
}

/**
 * Set the active flag of the bus
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @param flag flag for setting bus as active or inactive
 * @return false if no bus exists for idx
 */
extern "C" bool network_set_active_bus(networkWrapper *wnetwork, int idx, bool flag)
{
  idx--;
  return wnetwork->network->setActiveBus(idx,flag);
}

/**
 * Set the active flag of the branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param flag flag for setting bus as active or inactive
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_active_branch(networkWrapper *wnetwork, int idx, bool flag)
{
  idx--;
  return wnetwork->network->setActiveBranch(idx,flag);
}

/**
 * Clear the list of neighbors for the bus at idx
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @return false if no bus exists for idx
 */
extern "C" bool network_clear_branch_neighbors(networkWrapper *wnetwork, int idx)
{
  idx--;
  return wnetwork->network->clearBranchNeighbors(idx);
}

/**
 * Add local index for a branch attached to bus at idx
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @param br_idx local index of branch attached to bus
 * @return false if no bus exists for idx
 */
extern "C" bool network_add_branch_neighbor(networkWrapper *wnetwork, int idx, int br_idx)
{
  idx--;
  br_idx--;
  return wnetwork->network->addBranchNeighbor(idx,br_idx);
}

/**
 * Get status of the bus (local or ghosted)
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @return true if bus is locally held, false if it is ghosted 
 */
extern "C" bool network_get_active_bus(networkWrapper *wnetwork, int idx)
{
  idx--;
  return wnetwork->network->getActiveBus(idx);
}

/**
 * Get original index of the bus
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @return original index of bus 
 */
extern "C" int network_get_original_bus_index(networkWrapper *wnetwork, int idx)
{
  idx--;
  return wnetwork->network->getOriginalBusIndex(idx);
}

/**
 * Get global index of the bus
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @return global index of bus 
 */
extern "C" int network_get_global_bus_index(networkWrapper *wnetwork, int idx)
{
  idx--;
  int ret = wnetwork->network->getGlobalBusIndex(idx);
  ret++;
  return ret;
}

/**
 * Get status of the branch (local or ghosted)
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @return true if branch is locally held, false if it is ghosted 
 */
extern "C" bool network_get_active_branch(networkWrapper *wnetwork, int idx)
{
  idx--;
  return wnetwork->network->getActiveBranch(idx);
}

/**
 * Get global index of the branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @return global index of branch 
 */
extern "C" int network_get_global_branch_index(networkWrapper *wnetwork, int idx)
{
  idx--;
  int ret = wnetwork->network->getGlobalBranchIndex(idx);
  ret++;
  return ret;
}

/**
 * Get original indices of the two buses at each end of a branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param idx1 original index of "from" bus 
 * @param idx1 original index of "to" bus 
 */
extern "C" void network_get_original_branch_endpoints(networkWrapper *wnetwork,
    int idx, int *idx1, int *idx2)
{
  idx--;
  wnetwork->network->getOriginalBranchEndpoints(idx, idx1,idx2);
}

/**
 * Return the number of branches connected to a bus. This can be used to
 * allocate arrays to hold the branch indices
 * @param network pointer to Fortran network object
 * @param idx local bus index
 * @return the number of branches attached to this bus
 */
extern "C" int network_get_num_connected_branches(networkWrapper *wnetwork, int idx)
{
  idx--;
  return wnetwork->network->getConnectedBranches(idx).size();
}

/**
 * Return list of branches connected to bus
 * @param network pointer to Fortran network object
 * @param idx local bus index
 * @param branches array to hold list of branch indices
 */
extern "C" void network_get_connected_branches(networkWrapper *wnetwork, int idx, int *branches)
{
  int i;
  idx--;
  std::vector<int> t_branches = wnetwork->network->getConnectedBranches(idx);
  int size = t_branches.size();
  for (i=0; i<size; i++) {
    branches[i] = t_branches[i]+1;
  }
}

/**
 * Return list of buses connected to central bus via one branch
 * @param network pointer to Fortran network object
 * @param idx local bus index
 * @param buses array to hold list of bus indices
 */
extern "C" void network_get_connected_buses(networkWrapper *wnetwork, int idx, int *buses)
{
  int i;
  idx--;
  std::vector<int> t_buses = wnetwork->network->getConnectedBuses(idx);
  int size = t_buses.size();
  for (i=0; i<size; i++) {
    buses[i] = t_buses[i]+1;
  }
}

/**
 * Return indices of buses at either end of branch
 * @param network pointer to Fortran network object
 * @param idx local branch index
 * @param bus1 local index of bus at one end of branch
 * @param bus2 local index of bus at other end of branch
 */
extern "C" void network_get_branch_endpoints(networkWrapper *wnetwork, int idx,
    int *idx1, int *idx2)
{
  idx--;
  wnetwork->network->getBranchEndpoints(idx,idx1,idx2);
  *idx1 = *idx1+1;
  *idx2 = *idx2+1;
}

/**
 * Partition the network over the available processes
 * @param network pointer to Fortran network object
 */
extern "C" void network_partition(networkWrapper *wnetwork)
{
  wnetwork->network->partition();
}

/**
 * Clean all ghost buses and branches from the system. This can be used
 * before repartitioning the network. This operation also removes all exchange
 * buffers, so these need to be reallocated after calling this method
 * @param network pointer to Fortran network object
 */
extern "C" void network_clean(networkWrapper *wnetwork)
{
  wnetwork->network->clean();
}

/**
 * Allocate array of pointers to buffers for exchanging data for ghost buses
 * @param network GridPACK network object
 * @param isize size of exchange buffer (in bytes)
 */
extern "C" void network_alloc_xc_bus_pointers(networkWrapper *wnetwork, int
    isize)
{
  wnetwork->network->allocXCBusPointers(isize);
}

/**
 * Allocate array of pointers to buffers for exchanging data for ghost branches
 * @param network GridPACK network object
 * @param isize size of exchange buffer (in bytes)
 */
extern "C" void network_alloc_xc_branch_pointers(networkWrapper *wnetwork,
    int isize)
{
  wnetwork->network->allocXCBranchPointers(isize);
}

/**
 * Store location of externally allocated bus buffer within a network
 * @param p_network GridPACK network object
 * @param idx local index of bus associated with buffer
 * @param ptr location of buffer
 */
extern "C" void network_set_xc_bus_buffer(networkWrapper *wnetwork, int idx,
    void *data)
{
  idx--;
  wnetwork->network->setXCBusBuffer(idx, data);
}

/**
 * Store location of externally allocated branch buffer within a network
 * @param p_network GridPACK network object
 * @param idx local index of bus associated with buffer
 * @param ptr location of buffer
 */
extern "C" void network_set_xc_branch_buffer(networkWrapper *wnetwork,
    int idx, void *data)
{
  idx--;
  wnetwork->network->setXCBranchBuffer(idx, data);
}

/**
 * This function must be called before calling the update bus routine.
 * It initializes data structures for the bus update
 * @param network pointer to Fortran network object
 */
extern "C" void network_init_bus_update(networkWrapper *wnetwork)
{
  wnetwork->network->initBusUpdate();
}

/**
 * Update the bus ghost values. This is a
 * collective operation across all processors.
 * @param network pointer to Fortran network object
 */
extern "C" void network_update_buses(networkWrapper *wnetwork)
{
  wnetwork->network->updateBuses();
}

/**
 * This function must be called before calling the update branch routine.
 * It initializes data structures for the branch update
 * @param network pointer to Fortran network object
 */
extern "C" void network_init_branch_update(networkWrapper *wnetwork)
{
  wnetwork->network->initBranchUpdate();
}

/**
 * Update the branch ghost values. This is a
 * collective operation across all processors.
 * @param network pointer to Fortran network object
 */
extern "C" void network_update_branches(networkWrapper *wnetwork)
{
  wnetwork->network->updateBranches();
}

/**
 * Return a void pointer to the internal fortran bus object
 * @param network pointer to Fortran network object
 * @param idx index of bus object
 * @return void pointer to fortran bus object
 */
extern "C" void* network_get_bus(networkWrapper *wnetwork, int idx)
{
  idx--;
  return wnetwork->network->getBus(idx)->getFortranPointer();
}

/**
 * Return a void pointer to the internal fortran branch object
 * @param network pointer to Fortran network object
 * @param idx index of branch object
 * @return void pointer to fortran branch object
 */
extern "C" void* network_get_branch(networkWrapper *wnetwork, int idx)
{
  idx--;
  return wnetwork->network->getBranch(idx)->getFortranPointer();
}

/**
 * Return a void pointer to the internal fortran bus data collection object
 * @param network pointer to Fortran network object
 * @param idx index of bus object
 * @return void pointer to fortran bus data collection object
 */
extern "C" void* network_get_bus_data(networkWrapper *wnetwork, int idx)
{
  idx--;
  return wnetwork->network->getBusData(idx).get();
}

/**
 * Return a void pointer to the internal fortran branch data collection object
 * @param network pointer to Fortran network object
 * @param idx index of branch object
 * @return void pointer to fortran branch data collection object
 */
extern "C" void* network_get_branch_data(networkWrapper *wnetwork, int idx)
{
  idx--;
  return wnetwork->network->getBranchData(idx).get();
}
