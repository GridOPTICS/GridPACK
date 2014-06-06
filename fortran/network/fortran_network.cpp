#include "gridpack/parallel/communicator.hpp"
#include "fortran_network.hpp"

/**
 * Create a new network and pass an integer handle to it back to the calling
 * program
 */
extern "C" int p_create_network(void *comm_id)
{
  gridpack::parallel::Communicator *comm
    = static_cast<gridpack::parallel::Communicator*>(comm_id);
  int size = p_networks.size();
  // Find a p_network struct that isn't being used and initialize a new network
  int i, idx;
  i = 0;
  while (i<size) {
    if (!p_networks[i].active) {
      idx = i;
      p_networks[i].active = true;
      p_networks[i].network.reset(new FortranNetwork(*comm));
      break;
    }
    i++;
    if (i==size) {
      // TODO: Some kind of error
    }
  }
  return idx;
}

/**
 * Utility function to check that network handle is okay
 * @param n_handle network handle
 */
void p_checkHandle(int n_handle)
{
  if (n_handle < 0 || n_handle > MAX_NETWORKS) {
    printf("Handle out of bounds. handle: %d maximum networks: %d\n",
        n_handle, MAX_NETWORKS);
    // TODO: Some kind of error
  } else if (!p_networks[n_handle].active) {
    printf("Handle not active\n");
    // TODO: Some kind of error
  }
}


/**
 * Add a bus locally to the network
 * @param n_handle network handle
 * @param idx original index of bus
 */
extern "C" void p_add_bus(int n_handle, int idx)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->addBus(idx);
}

/**
 * Add a branch locally to the network. A branch is defined by
 * buses at either end
 * @param n_handle network handle
 * @param idx1 original bus index of bus 1
 * @param idx2 original bus index of bus 2
 */
extern "C" void p_add_branch(int n_handle, int idx1, int idx2)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->addBranch(idx1,idx2);
}

/**
 * Number of local buses (both active and inactive) on processor
 * @param n_handle network handle
 * @return number of buses
 */
extern "C" int p_num_buses(int n_handle)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->numBuses();
}

/**
 * Return the total number of buses in the entire network
 * @param n_handle network handle
 * @return total number of buses
 */
extern "C" int p_total_buses(int n_handle)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->totalBuses();
}

/**
 * Number of local branches (both active and inactive) on processor
 * @param n_handle network handle
 * @return number of branches
 */
extern "C" int p_num_branches(int n_handle)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->numBranches();
}

/**
 * Return the total number of branches in the entire network
 * @param n_handle network handle
 * @return total number of branches
 */
extern "C" int p_total_branches(int n_handle)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->totalBranches();
}

/**
 * Designate a bus as a reference bus.
 * @param n_handle network handle
 * @param idx local index of bus
 */
extern "C" void p_set_reference_bus(int n_handle, int idx)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->setReferenceBus(idx);
}

/**
 * Return index of reference bus.
 * @param n_handle network handle
 * @return local index of reference bus. If reference bus is not on this
 * processor then return -1.
 */
extern "C" void p_get_reference_bus(int n_handle)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->getReferenceBus();
}

/**
 * Set the original index of the bus (from configuration file)
 * @param n_handle network handle
 * @param idx local index of bus
 * @param o_idx original index assigned to bus
 * @return false if no bus exists for idx
 */
extern "C" bool p_set_original_bus_index(int n_handle, int idx, int o_idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->setOriginalBusIndex(idx,o_idx);
}

/**
 * Set the global index of the bus
 * @param n_handle network handle
 * @param idx local index of bus
 * @param g_idx global index to be assigned to bus
 * @return false if no bus exists for idx
 */
extern "C" bool p_set_global_bus_index(int n_handle, int idx, int g_idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->setGlobalBusIndex(idx,g_idx);
}

/**
 * Set the global index of the branch
 * @param n_handle network handle
 * @param idx local index of branch
 * @param g_idx global index to be assigned to branch
 * @return false if no branch exists for idx
 */
extern "C" bool p_set_global_branch_index(int n_handle, int idx, int g_idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->setGlobalBranchIndex(idx,g_idx);
}

/**
 * Set the original index of the bus at the "from" end of branch
 * @param n_handle network handle
 * @param idx local index of branch
 * @param b_idx original index of "from" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool p_set_original_bus_index1(int n_handle, int idx, int b_idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->setOriginalBusIndex1(idx,b_idx);
}

/**
 * Set the original index of the bus at the "to" end of branch
 * @param n_handle network handle
 * @param idx local index of branch
 * @param b_idx original index of "to" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool p_set_original_bus_index2(int n_handle, int idx, int b_idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->setOriginalBusIndex2(idx,b_idx);
}

/**
 * Set the global index of the bus at the "from" end of branch
 * @param n_handle network handle
 * @param idx local index of branch
 * @param b_idx global index of "from" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool p_set_global_bus_index1(int n_handle, int idx, int b_idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->setGlobalBusIndex1(idx,b_idx);
}

/**
 * Set the global index of the bus at the "to" end of branch
 * @param n_handle network handle
 * @param idx local index of branch
 * @param b_idx global index of "to" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool p_set_global_bus_index2(int n_handle, int idx, int b_idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->setGlobalBusIndex2(idx,b_idx);
}

/**
 * Set the local index of the bus at the "from" end of branch
 * @param n_handle network handle
 * @param idx local index of branch
 * @param b_idx local index of "from" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool p_set_local_bus_index1(int n_handle, int idx, int b_idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->setLocalBusIndex1(idx,b_idx);
}

/**
 * Set the local index of the bus at the "to" end of branch
 * @param n_handle network handle
 * @param idx local index of branch
 * @param b_idx local index of "to" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool p_set_local_bus_index2(int n_handle, int idx, int b_idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->setLocalBusIndex2(idx,b_idx);
}

/**
 * Set the active flag of the bus
 * @param n_handle network handle
 * @param idx local index of bus
 * @param flag flag for setting bus as active or inactive
 * @return false if no bus exists for idx
 */
extern "C" bool p_set_active_bus(int n_handle, int idx, bool flag)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->setActiveBus(idx,flag);
}

/**
 * Set the active flag of the branch
 * @param n_handle network handle
 * @param idx local index of branch
 * @param flag flag for setting bus as active or inactive
 * @return false if no branch exists for idx
 */
extern "C" bool p_set_active_branch(int n_handle, int idx, bool flag)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->setActiveBranch(idx,flag);
}

/**
 * Clear the list of neighbors for the bus at idx
 * @param n_handle network handle
 * @param idx local index of bus
 * @return false if no bus exists for idx
 */
extern "C" bool p_clear_branch_neighbors(int n_handle, int idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->clearBranchNeighbors(idx);
}

/**
 * Add local index for a branch attached to bus at idx
 * @param n_handle network handle
 * @param idx local index of bus
 * @param br_idx local index of branch attached to bus
 * @return false if no bus exists for idx
 */
extern "C" bool p_add_branch_neighbor(int n_handle, int idx, int br_idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->addBranchNeighbor(idx,br_idx);
}

/**
 * Get status of the bus (local or ghosted)
 * @param n_handle network handle
 * @param idx local index of bus
 * @return true if bus is locally held, false if it is ghosted 
 */
extern "C" bool p_get_active_bus(int n_handle, int idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->getActiveBus(idx);
}

/**
 * Get original index of the bus
 * @param n_handle network handle
 * @param idx local index of bus
 * @return original index of bus 
 */
extern "C" int p_get_original_bus_index(int n_handle, int idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->getOriginalBusIndex(idx);
}

/**
 * Get global index of the bus
 * @param n_handle network handle
 * @param idx local index of bus
 * @return global index of bus 
 */
extern "C" int p_get_global_bus_index(int n_handle, int idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->getGlobalBusIndex(idx);
}

/**
 * Get status of the branch (local or ghosted)
 * @param n_handle network handle
 * @param idx local index of branch
 * @return true if branch is locally held, false if it is ghosted 
 */
extern "C" bool p_get_active_branch(int n_handle, int idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->getActiveBranch(idx);
}

/**
 * Get global index of the branch
 * @param n_handle network handle
 * @param idx local index of branch
 * @return global index of branch 
 */
extern "C" int p_get_global_branch_index(int n_handle, int idx)
{
  p_checkHandle(n_handle);
  return p_networks[n_handle].network->getGlobalBranchIndex(idx);
}

/**
 * Get original indices of the two buses at each end of a branch
 * @param n_handle network handle
 * @param idx local index of branch
 * @param idx1 original index of "from" bus 
 * @param idx1 original index of "to" bus 
 */
extern "C" void p_get_original_branch_endpoints(int n_handle, int idx,
                                                int *idx1, int *idx2)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->getOriginalBranchEndpoints(idx,
      idx1,idx2);
}

/**
 * Return the number of branches connected to a bus. This can be used to
 * allocate arrays to hold the branch indices
 * @param n_handle network handle
 * @param idx local bus index
 * @return the number of branches attached to this bus
 */
extern "C" int p_num_connected_branches(int n_handle, int idx)
{
  p_checkHandle(n_handle);
  if (idx<0 || idx >= p_networks[n_handle].network->numBuses()) {
    // TODO: Some kind of error
    return -1;
  } else {
    return p_networks[n_handle].network->getConnectedBranches(idx).size();
  }
}

/**
 * Return list of branches connected to bus
 * @param n_handle network handle
 * @param idx local bus index
 * @param branches array to hold list of branch indices
 */
extern "C" void p_get_connected_branches(int n_handle, int idx, int *branches)
{
  p_checkHandle(n_handle);
  int i;
  std::vector<int> t_branches =
    p_networks[n_handle].network->getConnectedBranches(idx);
  int size = t_branches.size();
  for (i=0; i<size; i++) {
    branches[i] = t_branches[i];
  }
}

/**
 * Return list of buses connected to central bus via one branch
 * @param n_handle network handle
 * @param idx local bus index
 * @param buses array to hold list of bus indices
 */
extern "C" void p_get_connected_buses(int n_handle, int idx, int *buses)
{
  p_checkHandle(n_handle);
  if (idx<0 || idx >= p_networks[n_handle].network->numBuses()) {
    // TODO: Some kind of error
  } else {
    int i;
    std::vector<int> t_buses =
      p_networks[n_handle].network->getConnectedBuses(idx);
    int size = t_buses.size();
    for (i=0; i<size; i++) {
      buses[i] = t_buses[i];
    }
  }
}

/**
 * Return indices of buses at either end of branch
 * @param n_handle network handle
 * @param idx local branch index
 * @param bus1 local index of bus at one end of branch
 * @param bus2 local index of bus at other end of branch
 */
extern "C" void p_get_branch_endpoints(int n_handle, int idx, int *idx1, int *idx2)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->getBranchEndpoints(idx,idx1,idx2);
}

/**
 * Partition the network over the available processes
 * @param n_handle network handle
 */
extern "C" void p_partition(int n_handle)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->partition();
}

/**
 * Clean all ghost buses and branches from the system. This can be used
 * before repartitioning the network. This operation also removes all exchange
 * buffers, so these need to be reallocated after calling this method
 * @param n_handle network handle
 */
extern "C" void p_clean(int n_handle)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->clean();
}

/**
 * Allocate buffers for exchanging data for ghost buses
 * @param size size of buffers that will be assigned to pointers
 * @param n_handle network handle
 */
extern "C" void p_alloc_xc_bus_pointers(int n_handle, int size)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->allocXCBus(size);
}

/**
 * Store location of externally allocated buffer within a network
 * @param n_handle network handle
 * @param idx local index of bus associated with buffer
 * @param ptr location of buffer
 */
extern "C" void p_set_xc_bus_buffer(int n_handle, int idx, void *ptr)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->setXCBusBuffer(idx, ptr);
}

/**
 * Allocate buffers for exchanging data for ghost branches
 * @param size size of buffers that will be assigned to pointers
 * @param n_handle network handle
 */
extern "C" void p_alloc_xc_branch_pointers(int n_handle, int size)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->allocXCBranch(size);
}

/**
 * Store location of externally allocated buffer within a network
 * @param n_handle network handle
 * @param idx local index of branch associated with buffer
 * @param ptr location of buffer
 */
extern "C" void p_set_xc_branch_buffer(int n_handle, int idx, void *ptr)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->setXCBranchBuffer(idx, ptr);
}

/**
 * This function must be called before calling the update bus routine.
 * It initializes data structures for the bus update
 * @param n_handle network handle
 */
extern "C" void p_init_bus_update(int n_handle)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->initBusUpdate();
}

/**
 * Update the bus ghost values. This is a
 * collective operation across all processors.
 * @param n_handle network handle
 */
extern "C" void p_update_buses(int n_handle)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->updateBuses();
}

/**
 * This function must be called before calling the update branch routine.
 * It initializes data structures for the branch update
 * @param n_handle network handle
 */
extern "C" void p_init_branch_update(int n_handle)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->initBranchUpdate();
}

/**
 * Update the branch ghost values. This is a
 * collective operation across all processors.
 * @param n_handle network handle
 */
extern "C" void p_update_branches(int n_handle)
{
  p_checkHandle(n_handle);
  p_networks[n_handle].network->updateBranches();
}
