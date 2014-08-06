#include "gridpack/parallel/communicator.hpp"
#include "gridpack/network/base_network.hpp"
#include "../component/fortran_component.hpp"
//#include "fortran_network.hpp"

typedef gridpack::network::BaseNetwork<
    gridpack::fortran_component::FortranBusComponent,
    gridpack::fortran_component::FortranBranchComponent>
    FortranNetwork;

/**
 * Create a new network
 * @param network pointer to new Fortran network object
 * @param comm communicator
 */
extern "C" void network_create(FortranNetwork **network,
    gridpack::parallel::Communicator *comm)
{
  *network = new FortranNetwork(*comm);
}

/**
 * Destroy old network
 * @param network pointer to Fortran network object
 */
extern "C" void network_destroy(FortranNetwork **network)
{
  delete (*network);
  *network = NULL;
}
#if 0
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
#endif


/**
 * Add a bus locally to the network
 * @param network pointer to Fortran network object
 * @param idx original index of bus
 */
extern "C" void network_add_bus(FortranNetwork *network, int idx)
{
  network->addBus(idx);
}

/**
 * Add a branch locally to the network. A branch is defined by
 * buses at either end
 * @param network pointer to Fortran network object
 * @param idx1 original bus index of "from" bus
 * @param idx2 original bus index of "to" bus
 */
extern "C" void network_add_branch(FortranNetwork *network, int idx1, int idx2)
{
  network->addBranch(idx1,idx2);
}

/**
 * Number of local buses (both active and inactive) on processor
 * @param network pointer to Fortran network object
 * @return number of buses
 */
extern "C" int network_num_buses(FortranNetwork *network)
{
  return network->numBuses();
}

/**
 * Return the total number of buses in the entire network
 * @param network pointer to Fortran network object
 * @return total number of buses
 */
extern "C" int network_total_buses(FortranNetwork *network)
{
  return network->totalBuses();
}

/**
 * Number of local branches (both active and inactive) on processor
 * @param network pointer to Fortran network object
 * @return number of branches
 */
extern "C" int network_num_branches(FortranNetwork *network)
{
  return network->numBranches();
}

/**
 * Return the total number of branches in the entire network
 * @param network pointer to Fortran network object
 * @return total number of branches
 */
extern "C" int network_total_branches(FortranNetwork *network)
{
  return network->totalBranches();
}

/**
 * Designate a bus as a reference bus.
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 */
extern "C" void network_set_reference_bus(FortranNetwork *network, int idx)
{
  network->setReferenceBus(idx);
}

/**
 * Return index of reference bus.
 * @param network pointer to Fortran network object
 * @return local index of reference bus. If reference bus is not on this
 * processor then return -1.
 */
extern "C" int network_get_reference_bus(FortranNetwork *network)
{
  return network->getReferenceBus();
}

/**
 * Set the original index of the bus (from configuration file)
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @param o_idx original index assigned to bus
 * @return false if no bus exists for idx
 */
extern "C" bool network_set_original_bus_index(FortranNetwork *network, int idx, int o_idx)
{
  return network->setOriginalBusIndex(idx,o_idx);
}

/**
 * Set the global index of the bus
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @param g_idx global index to be assigned to bus
 * @return false if no bus exists for idx
 */
extern "C" bool network_set_global_bus_index(FortranNetwork *network, int idx, int g_idx)
{
  return network->setGlobalBusIndex(idx,g_idx);
}

/**
 * Set the global index of the branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param g_idx global index to be assigned to branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_global_branch_index(FortranNetwork *network, int idx, int g_idx)
{
  return network->setGlobalBranchIndex(idx,g_idx);
}

/**
 * Set the original index of the bus at the "from" end of branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param b_idx original index of "from" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_original_bus_index1(FortranNetwork *network, int idx, int b_idx)
{
  return network->setOriginalBusIndex1(idx,b_idx);
}

/**
 * Set the original index of the bus at the "to" end of branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param b_idx original index of "to" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_original_bus_index2(FortranNetwork *network, int idx, int b_idx)
{
  return network->setOriginalBusIndex2(idx,b_idx);
}

/**
 * Set the global index of the bus at the "from" end of branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param b_idx global index of "from" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_global_bus_index1(FortranNetwork *network, int idx, int b_idx)
{
  return network->setGlobalBusIndex1(idx,b_idx);
}

/**
 * Set the global index of the bus at the "to" end of branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param b_idx global index of "to" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_global_bus_index2(FortranNetwork *network, int idx, int b_idx)
{
  return network->setGlobalBusIndex2(idx,b_idx);
}

/**
 * Set the local index of the bus at the "from" end of branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param b_idx local index of "from" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_local_bus_index1(FortranNetwork *network, int idx, int b_idx)
{
  return network->setLocalBusIndex1(idx,b_idx);
}

/**
 * Set the local index of the bus at the "to" end of branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param b_idx local index of "to" bus for this branch
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_local_bus_index2(FortranNetwork *network, int idx, int b_idx)
{
  return network->setLocalBusIndex2(idx,b_idx);
}

/**
 * Set the active flag of the bus
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @param flag flag for setting bus as active or inactive
 * @return false if no bus exists for idx
 */
extern "C" bool network_set_active_bus(FortranNetwork *network, int idx, bool flag)
{
  return network->setActiveBus(idx,flag);
}

/**
 * Set the active flag of the branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param flag flag for setting bus as active or inactive
 * @return false if no branch exists for idx
 */
extern "C" bool network_set_active_branch(FortranNetwork *network, int idx, bool flag)
{
  return network->setActiveBranch(idx,flag);
}

/**
 * Clear the list of neighbors for the bus at idx
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @return false if no bus exists for idx
 */
extern "C" bool network_clear_branch_neighbors(FortranNetwork *network, int idx)
{
  return network->clearBranchNeighbors(idx);
}

/**
 * Add local index for a branch attached to bus at idx
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @param br_idx local index of branch attached to bus
 * @return false if no bus exists for idx
 */
extern "C" bool network_add_branch_neighbor(FortranNetwork *network, int idx, int br_idx)
{
  return network->addBranchNeighbor(idx,br_idx);
}

/**
 * Get status of the bus (local or ghosted)
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @return true if bus is locally held, false if it is ghosted 
 */
extern "C" bool network_get_active_bus(FortranNetwork *network, int idx)
{
  return network->getActiveBus(idx);
}

/**
 * Get original index of the bus
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @return original index of bus 
 */
extern "C" int network_get_original_bus_index(FortranNetwork *network, int idx)
{
  return network->getOriginalBusIndex(idx);
}

/**
 * Get global index of the bus
 * @param network pointer to Fortran network object
 * @param idx local index of bus
 * @return global index of bus 
 */
extern "C" int network_get_global_bus_index(FortranNetwork *network, int idx)
{
  return network->getGlobalBusIndex(idx);
}

/**
 * Get status of the branch (local or ghosted)
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @return true if branch is locally held, false if it is ghosted 
 */
extern "C" bool network_get_active_branch(FortranNetwork *network, int idx)
{
  return network->getActiveBranch(idx);
}

/**
 * Get global index of the branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @return global index of branch 
 */
extern "C" int network_get_global_branch_index(FortranNetwork *network, int idx)
{
  return network->getGlobalBranchIndex(idx);
}

/**
 * Get original indices of the two buses at each end of a branch
 * @param network pointer to Fortran network object
 * @param idx local index of branch
 * @param idx1 original index of "from" bus 
 * @param idx1 original index of "to" bus 
 */
extern "C" void network_get_original_branch_endpoints(FortranNetwork *network,
    int idx, int *idx1, int *idx2)
{
  network->getOriginalBranchEndpoints(idx, idx1,idx2);
}

/**
 * Return the number of branches connected to a bus. This can be used to
 * allocate arrays to hold the branch indices
 * @param network pointer to Fortran network object
 * @param idx local bus index
 * @return the number of branches attached to this bus
 */
extern "C" int network_get_num_connected_branches(FortranNetwork *network, int idx)
{
  return network->getConnectedBranches(idx).size();
}

/**
 * Return list of branches connected to bus
 * @param network pointer to Fortran network object
 * @param idx local bus index
 * @param branches array to hold list of branch indices
 */
extern "C" void network_get_connected_branches(FortranNetwork *network, int idx, int *branches)
{
  int i;
  std::vector<int> t_branches = network->getConnectedBranches(idx);
  int size = t_branches.size();
  for (i=0; i<size; i++) {
    branches[i] = t_branches[i];
  }
}

/**
 * Return list of buses connected to central bus via one branch
 * @param network pointer to Fortran network object
 * @param idx local bus index
 * @param buses array to hold list of bus indices
 */
extern "C" void network_get_connected_buses(FortranNetwork *network, int idx, int *buses)
{
  int i;
  std::vector<int> t_buses = network->getConnectedBuses(idx);
  int size = t_buses.size();
  for (i=0; i<size; i++) {
    buses[i] = t_buses[i];
  }
}

/**
 * Return indices of buses at either end of branch
 * @param network pointer to Fortran network object
 * @param idx local branch index
 * @param bus1 local index of bus at one end of branch
 * @param bus2 local index of bus at other end of branch
 */
extern "C" void network_get_branch_endpoints(FortranNetwork *network, int idx, int *idx1, int *idx2)
{
  network->getBranchEndpoints(idx,idx1,idx2);
}

/**
 * Partition the network over the available processes
 * @param network pointer to Fortran network object
 */
extern "C" void network_partition(FortranNetwork *network)
{
  network->partition();
}

/**
 * Clean all ghost buses and branches from the system. This can be used
 * before repartitioning the network. This operation also removes all exchange
 * buffers, so these need to be reallocated after calling this method
 * @param network pointer to Fortran network object
 */
extern "C" void network_clean(FortranNetwork *network)
{
  network->clean();
}

#if 0
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
#endif

/**
 * This function must be called before calling the update bus routine.
 * It initializes data structures for the bus update
 * @param n_handle network handle
 */
extern "C" void network_init_bus_update(FortranNetwork *network)
{
  network->initBusUpdate();
}

/**
 * Update the bus ghost values. This is a
 * collective operation across all processors.
 * @param n_handle network handle
 */
extern "C" void network_update_buses(FortranNetwork *network)
{
  network->updateBuses();
}

/**
 * This function must be called before calling the update branch routine.
 * It initializes data structures for the branch update
 * @param n_handle network handle
 */
extern "C" void network_init_branch_update(FortranNetwork *network)
{
  network->initBranchUpdate();
}

/**
 * Update the branch ghost values. This is a
 * collective operation across all processors.
 * @param n_handle network handle
 */
extern "C" void network_update_branches(FortranNetwork *network)
{
  network->updateBranches();
}
