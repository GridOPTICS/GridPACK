// -------------------------------------------------------------
/**
 * @file   base_network.hpp
 * @author Bruce Palmer, William Perkins
 * @date   April 3, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _base_network_h_
#define _base_network_h_

#include <vector>
#include <map>
#include <boost/smart_ptr/shared_ptr.hpp>
#include "gridpack/parallel/distributed.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"

// -------------------------------------------------------------
//  class BaseNetwork:
//  This is the base class for creating distributed networks. It
//  is basically something that supports the network topology,
//  allows fields to be added and subtracted to the buses (nodes)
//  and branches (edges) of the network, and implements ghost
//  updates on the network. The BaseNetwork class does not
//  contain the partitioner but it does contain methods that
//  allow the partitioner to move grid elements around and
//  create ghost elements.
// -------------------------------------------------------------
namespace gridpack {
namespace network {

template <class _bus, class _branch>
class BaseNetwork {
public:
/**
 * Default constructor.
 */
BaseNetwork(void)
{
  //TODO: Get default parallel configuration for world group
  p_refBus = -1;
}

/**
 * Constructor with parallel configuration for running on
 * subgroups or communicators other than the world
 * communicator
 */
#if 0
BaseNetwork::BaseNetwork(ParallelEnv configuration)
{
  p_configuration = configuration;
  p_refBus = -1;
}
#endif

/**
 * Default destructor.
 */
~BaseNetwork(void)
{
}

/**
 * Add a bus locally to the network
 * @param idx: original index of  bus
 */
void addBus(int idx)
{
  p_originalBusIndex.push_back(idx);
  p_globalBusIndex.push_back(-1);
  p_activeBus.push_back(true);
  boost::shared_ptr<_bus> (new _bus) bus;
  p_buses.push_back(bus);
  boost::shared_ptr<gridpack::component::DataCollection>(new
      gridpack::component::DataCollection) data;
  p_busData.push_back(data);
}

/**
 * Add a branch locally to the network. A branch is defined by
 * buses at either end
 * @param idx1: original bus index of bus 1
 * @param idx2: original bus index of bus 2
 */
void addBranch(int idx1, int idx2)
{
  p_originalBranchIndex1.push_back(idx1);
  p_originalBranchIndex2.push_back(idx2);
  p_globalBranchIndex1.push_back(-1);
  p_globalBranchIndex2.push_back(-1);
  p_activeBranch.push_back(true);
  boost::shared_ptr<_branch> (new _branch) branch;
  p_branches.push_back(branch);
  boost::shared_ptr<gridpack::component::DataCollection>(new
      gridpack::component::DataCollection) data;
  p_branchData.push_back(data);
}

/**
 * Designate a bus as a reference bus.
 * @param idx: local index of bus
 */
void setReferenceBus(int idx)
{
  p_refBus = idx;
}

/**
 * Return index of reference bus.
 * @return: local index of reference bus. If reference bus is not on this
 * processor then return -1.
 */
int getReferenceBus(void) const
{
  return p_refBus;
}

/**
 * Set the global index of the bus
 * @param idx: local index of bus
 * @param g_idx: global index to be assigned to bus
 * @return: false if no bus exists for idx
 */
bool setGlobalBusIndex(int idx, int g_idx)
{
  if (idx < 0 || idx >= p_globalBusIndex.size()) {
    return false;
  } else {
    p_globalBusIndex[idx] = g_idx;
    return true;
  }
}

/**
 * Set the global index of the branch
 * @param idx: local index of branch
 * @param g_idx: global index to be assigned to branch
 * @return: false if no branch exists for idx
 */
bool setGlobalBranchIndex(int idx, int g_idx)
{
  if (idx < 0 || idx >= p_globalBranchIndex.size()) {
    return false;
  } else {
    p_globalBranchIndex[idx] = g_idx;
    return true;
  }
}

/**
 * Set the active flag of the bus
 * @param idx: local index of bus
 * @param flag: flag for setting bus as active or inactive
 * @return: false if no bus exists for idx
 */
bool setActiveBus(int idx, bool flag)
{
  if (idx < 0 || idx >= p_activeBus.size()) {
    return false;
  } else {
    p_activeBus[idx] = flag;
    return true;
  }
}

/**
 * Set the active flag of the branch
 * @param idx: local index of branch
 * @param flag: flag for setting bus as active or inactive
 * @return: false if no branch exists for idx
 */
bool setActiveBranch(int idx, bool flag)
{
  if (idx < 0 || idx >= p_activeBranch.size()) {
    return false;
  } else {
    p_activeBranch[idx] = flag;
    return true;
  }
}

/**
 * Retrieve a pointer to an existing bus
 * @param idx: local index of requested bus
 * @return: a pointer to the requested bus.
 */
boost::shared_ptr<_bus> getBus(int idx)
{
  if (idx<0 || idx >= p_buses->size()) {
    // TODO: Some kind of error
  } else {
    return p_buses[idx];
  }
}

/**
 * Retrieve a pointer to an existing branch
 * @param idx: local index of requested branch
 * @return: a pointer to the requested branch
 */
boost::shared_ptr<_branch> getBranch(int idx)
{
  if (idx<0 || idx >= p_branches->size()) {
    // TODO: Some kind of error
  } else {
    return p_branches[idx];
  }
}

/**
 * Retrieve a pointer to the DataCollection object associated with bus indexed
 * by idx
 * @param idx: local index of requested bus
 * @return: a pointer to the requested bus data
 */
boost::shared_ptr<gridpack::component::DataCollection> getBusData(int idx)
{
  if (idx<0 || idx >= p_buses->size()) {
    // TODO: Some kind of error
  } else {
    return p_busData[idx];
  }
}

/**
 * Retrieve a pointer to the DataCollection object associated with branch
 * indexed by idx
 * @param idx: local index of requested branch
 * @return: a pointer to the requested branch data
 */
boost::shared_ptr<gridpack::component::DataCollection> getBranchData(int idx)
{
  if (idx<0 || idx >= p_branches->size()) {
    // TODO: Some kind of error
  } else {
    return p_branchData[idx];
  }
}

#if 0
/**
 * Delete an existing bus field
 * @param name: a string representing the name of the field
 *       to be deleted
 */
void deleteBusField(std::string name)
{
  std::map<std::string,
boost::shared_ptr<gridpack::network::BusField<gridpack::component::BaseNetworkComponent> > >::iterator bus;
  bus = p_busFields.find(name);
  if (bus != p_busFields.end()) {
    p_busFields.erase(bus);
  }
}

/**
 * Delete an existing branch field
 * @param name: a string representing the name of the field
 *       to be deleted
 */
void deleteBranchField(std::string name)
{
  std::map<std::string,
boost::shared_ptr<gridpack::network::BranchField<gridpack::component::BaseNetworkComponent> > >::iterator branch;
  branch = p_branchFields.find(name);
  if (branch != p_branchFields.end()) {
    p_branchFields.erase(branch);
  }
}
#endif

/**
 * Return list of branches connected to bus
 * @param bus: local bus index
 * @return: vector of local branch indices
 */
std::vector<int> getConnectedBranches(int idx) const
{
  return p_branchNeighbors[idx];
}

/**
 * Return list of buses connected to central bus via one branch
 * @param bus: local bus index
 * @return: vector of local bus indices
 */
std::vector<int> getConnectedBuses(int idx) const
{
  std::vector<int> branches = p_branchNeighbors[idx];
  int size = branches.size();
  std::vector<int> ret;
  int i, j;
  for (i=0; i<size; i++) {
    j = branches[i];
    if (p_localBranchIndex1[j] != idx) {
      ret.push_back(p_localBranchIndex1[j]);
    } else {
      ret.push_back(p_localBranchIndex2[j]);
    }
  }
  return ret;
}

/**
 * Return indices of buses at either end of branch
 * @param bus: local branch index
 * @param bus1: local index of bus at one end of branch
 * @param bus2: local index of bus at other end of branch
 */
void getBusEndpoints(int idx, int *bus1, int *bus2) const
{
  *bus1 = p_localBranchIndex1[idx];
  *bus2 = p_localBranchIndex2[idx];
}

/**
 * Clean all ghost buses and branches from the system. This can be used
 * before repartitioning the network
 */
void clean(void)
{
  std::map<int, int> buses;
  std::map<int, int> branches;
  std::map<int, int>::iterator p;
  int i, j;

  // remove inactive branches
  int size = p_activeBranch.size();
  int new_id = 0;
  for (i=0; i<size; i++) {
    if (p_activeBranch[i]) {
      p_activeBranch[new_id] = p_activeBranch[i];
      p_originalBranchIndex1[new_id] = p_originalBranchIndex1[i];
      p_originalBranchIndex2[new_id] = p_originalBranchIndex2[i];
      p_globalBranchIndex1[new_id] = p_globalBranchIndex1[i];
      p_globalBranchIndex2[new_id] = p_globalBranchIndex2[i];
      p_localBranchIndex1[new_id] = p_localBranchIndex1[i];
      p_localBranchIndex2[new_id] = p_localBranchIndex2[i];
      p_branches[new_id] = p_branches[i];
      branches.insert(std::pair<int, int>(i,new_id));
      new_id++;
    }
  }
  // clean up the ends of the branch vectors
  i = size - new_id;
  for (i=0; i<size; i++) {
    p_activeBranch.pop_back();
    p_originalBranchIndex1.pop_back();
    p_originalBranchIndex2.pop_back();
    p_globalBranchIndex1.pop_back();
    p_globalBranchIndex2.pop_back();
    p_localBranchIndex1.pop_back();
    p_localBranchIndex2.pop_back();
    p_branches.pop_back();
  }

  // remove inactive buses
  size = p_activeBus.size();
  new_id = 0;
  for (i=0; i<size; i++) {
    if (p_activeBus[i]) {
      p_activeBus[new_id] = p_activeBus[i];
      p_originalBusIndex[new_id] = p_originalBusIndex[i];
      p_branchNeighbors[new_id] = p_branchNeighbors[i];
      p_busNeighbors[new_id] = p_busNeighbors[i];
      p_buses[new_id] = p_buses[i];
      buses.insert(std::pair<int, int>(i,new_id));
      new_id++;
    }
  }
  // clean up the ends of the bus vectors
  i = size - new_id;
  for (i=0; i<size; i++) {
    p_activeBus.pop_back();
    p_originalBusIndex.pop_back();
    p_globalBusIndex.pop_back();
    p_branchNeighbors.pop_back();
    p_busNeighbors.pop_back();
    p_buses.pop_back();
  }

  // reset all local indices
  size = p_activeBranch.size();
  for (i=0; i<size; i++) {
    p = buses.find(p_localBranchIndex1[i]);
    if (p != buses.end()) {
      p_localBranchIndex1[i] = p->second;
    } else {
      p_localBranchIndex1[i] = -1;
    }
    p = buses.find(p_localBranchIndex2[i]);
    if (p != buses.end()) {
      p_localBranchIndex2[i] = p->second;
    } else {
      p_localBranchIndex2[i] = -1;
    }
  }

  int jsize;
  size = p_activeBus.size();
  for (i=0; i<size; i++) {
    std::vector<int> neighbors = p_branchNeighbors[i];
    jsize = neighbors.size();
    p_branchNeighbors[i].clear();
    for (j=0; j<jsize; j++) {
      p = branches.find(neighbors[j]);
      if (p != branches.end()) {
        p_branchNeighbors[i].push_back(p->second);
      }
    }
    std::vector<int> neighbors = p_busNeighbors[i];
    jsize = neighbors.size();
    p_busNeighbors[i].clear();
    for (j=0; j<jsize; j++) {
      p = buses.find(neighbors[j]);
      if (p != branches.end()) {
        p_busNeighbors[i].push_back(p->second);
      }
    }
  }
  if (p_refBus != -1) {
    p_refBus = buses[p_refBus];
  }
}

/**
 * Update the ghost values of this field. This is a
 * collective operation across all processors.
 * @param field: name of the network field that must be
 *       updated
 */
void updateField(char *field)
{
}

protected:

/**
 * Protected copy constructor to avoid unwanted copies.
 */
BaseNetwork(const BaseNetwork& old)
{
}

private:

  /**
   * Vector that distinguishes active (local) buses from inactive (ghost) buses
   */
  std::vector<bool> p_activeBus;

  /**
   * Vector that distinguishes active (local) branches from inactive (ghost)
   * branches
   */
  std::vector<bool> p_activeBranch;

  /**
   * Original index of a bus (from the network topology file)
   */
  std::vector<int> p_originalBusIndex;

  /**
   * Unique global index of a bus assigned by the partitioner
   */
  std::vector<int> p_globalBusIndex;

  /**
   * Local indices of branches that are connected to a local bus
   */
  std::vector<std::vector<int> > p_branchNeighbors;

  /**
   * Local indices of buses that are connected to a local bus via a single
   * branch
   */
  std::vector<std::vector<int> > p_busNeighbors;

  /**
   * Original index of bus at one end of a branch
   */
  std::vector<int> p_originalBranchIndex1;

  /**
   * Original index of bus at other end of a branch
   */
  std::vector<int> p_originalBranchIndex2;

  /**
   * Global index of bus at one end of a branch
   */
  std::vector<int> p_globalBranchIndex1;

  /**
   * Global index of bus at other end of a branch
   */
  std::vector<int> p_globalBranchIndex2;

  /**
   * Local index of bus at one end of a branch
   */
  std::vector<int> p_localBranchIndex1;

  /**
   * Local index of bus at one end of a branch
   */
  std::vector<int> p_localBranchIndex2;

  /**
   * vector of bus objects
   */
  std::vector<boost::shared_ptr<_bus> > p_buses;

  /**
   * vector of branch objects
   */
  std::vector<boost::shared_ptr<_branch> > p_branches;

  /**
   * Parallel environment for network
   */
#if 0
  ParallelEnv p_configuration;
#endif

  int p_refBus;
  std::vector<boost::shared_ptr<gridpack::component::DataCollection> > p_busData;
  std::vector<boost::shared_ptr<gridpack::component::DataCollection> > p_branchData;
};
}  //namespace network
}  //namespace gridpack

#endif
