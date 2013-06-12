// -------------------------------------------------------------
/**
 * @file   base_network.hpp
 * @author Bruce Palmer, William Perkins
 * @date   2013-06-04 12:42:01 d3g096
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
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"

namespace gridpack {
namespace network {

// -------------------------------------------------------------
// A simple data class to assemble all bus related elements in
// the network into a single struct.
// -------------------------------------------------------------
template <class _bus>
class BusData {
  public:

/**
 *  Default constructor
 */
BusData(void)
{
  p_bus = boost::shared_ptr<_bus> (new _bus);
  p_activeBus = true;
  p_refFlag = false;
}

/**
 *  Default destructor
 */
~BusData(void)
{
}

/**
 *  Assigment operator
 */
BusData<_bus> & operator=(const BusData<_bus> & rhs)
{
  if (this == &rhs) return *this;
  p_activeBus = rhs.p_activeBus;
  p_originalBusIndex = rhs.p_originalBusIndex;
  p_globalBusIndex = rhs.p_globalBusIndex;
  p_branchNeighbors = rhs.p_branchNeighbors;
  p_bus = rhs.p_bus;
  p_data = rhs.p_data;
  p_refFlag = rhs.p_refFlag;
}

/**
 * Data elements on bus
 * p_activeBus: flag to identify buses that are "owned" by processor
 * p_originalBusIndex: original index of a bus (from the network topology file)
 * p_globalBusIndex: unique global index of a bus assigned by the partitioner
 * p_branchNeighbors: local indices of branches that are connected to a local bus
 * p_bus: pointer to bus object
 * p_data: pointer to data collection object
 * p_refFlag: true if this bus is the reference bus
 */
  bool                                                   p_activeBus;
  int                                                    p_originalBusIndex;
  int                                                    p_globalBusIndex;
  std::vector<int>                                       p_branchNeighbors;
  boost::shared_ptr<_bus>                                p_bus;
  boost::shared_ptr<gridpack::component::DataCollection> p_data;
  bool                                                   p_refFlag;
};

// -------------------------------------------------------------
// A simple data class to assemble all branch related elements
// in the network into a single struct.
// -------------------------------------------------------------
template <class _branch>
class BranchData {
  public:

/**
 *  Default constructor
 */
BranchData(void)
{
  p_branch = boost::shared_ptr<_branch> (new _branch);
  p_activeBranch = true;
}

/**
 *  Default destructor
 */
~BranchData(void)
{
}

/**
 *  Assigment operator
 */
BranchData<_branch> & operator=(const BranchData<_branch> & rhs)
{
  if (this == &rhs) return *this;
  p_activeBranch = rhs.p_activeBranch;
  p_globalBranchIndex = rhs.p_globalBranchIndex;
  p_originalBusIndex1 = rhs.p_originalBusIndex1;
  p_originalBusIndex2 = rhs.p_originalBusIndex2;
  p_globalBusIndex1 = rhs.p_globalBusIndex1;
  p_globalBusIndex2 = rhs.p_globalBusIndex2;
  p_localBusIndex1 = rhs.p_localBusIndex1;
  p_localBusIndex2 = rhs.p_localBusIndex2;
  p_branch = rhs.p_branch;
  p_data = rhs.p_data;
}

/**
 * Data elements on branch
 * p_activeBranch: flag to identify buses that are "owned" by processor
 * p_globalBranchIndex: unique global identifier for this branch
 * p_originalBusIndex1: original index of bus at "from" end of branch
 *      (from the network topology file)
 * p_originalBusIndex2: original index of bus at "to" end of branch
 *      (from the network topology file)
 * p_globalBusIndex1: global index of bus at "from" end of branch
 * p_globalBusIndex2: global index of bus at "to" end of branch
 * p_localBusIndex1: local index of bus at "from" end of branch
 * p_localBusIndex2: local index of bus at "to" end of branch
 * p_branch: pointer to branch object
 * p_data: pointer to data collection object
 */
  bool                                                   p_activeBranch;
  int                                                    p_globalBranchIndex;
  int                                                    p_originalBusIndex1;
  int                                                    p_originalBusIndex2;
  int                                                    p_globalBusIndex1;
  int                                                    p_globalBusIndex2;
  int                                                    p_localBusIndex1;
  int                                                    p_localBusIndex2;
  boost::shared_ptr<_branch>                             p_branch;
  boost::shared_ptr<gridpack::component::DataCollection> p_data;
};

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
  boost::shared_ptr<gridpack::network::BusData<_bus> >
    bus(new gridpack::network::BusData<_bus>);
  bus->p_originalBusIndex = idx;
  bus->p_globalBusIndex = -1;
  p_buses.push_back(bus);
}

/**
 * Add a branch locally to the network. A branch is defined by
 * buses at either end
 * @param idx1: original bus index of bus 1
 * @param idx2: original bus index of bus 2
 */
void addBranch(int idx1, int idx2)
{
  boost::shared_ptr<gridpack::network::BranchData<_branch> >
    branch(new gridpack::network::BranchData<_branch>);
  branch->p_originalBusIndex1 = idx1;
  branch->p_originalBusIndex2 = idx2;
  branch->p_globalBusIndex1 = -1;
  branch->p_globalBusIndex2 = -1;
  p_branches.push_back(branch);
}

/**
 * Number of local buses (both active and inactive) on processor
 * @return: number of buses
 */
int numBuses(void)
{
  return p_buses.size();
}

/**
 * Number of local branches (both active and inactive) on processor
 * @return: number of branches
 */
int numBranches(void)
{
  return p_branches.size();
}

/**
 * Designate a bus as a reference bus.
 * @param idx: local index of bus
 */
void setReferenceBus(int idx)
{
  p_refBus = idx;
  if (idx > 0 && idx <= p_buses.size()) {
    p_buses[idx]->p_refFlag = true;
  } else {
    p_buses[idx]->p_refFlag = false;
  }
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

// Bus and Branch modifiers

/**
 * Set the original index of the bus (from configuration file)
 * @param idx: local index of bus
 * @param o_idx: original index assigned to bus
 * @return: false if no bus exists for idx
 */
bool setOriginalBusIndex(int idx, int o_idx)
{
  if (idx < 0 || idx >= p_buses.size()) {
    return false;
  } else {
    p_buses[idx]->p_originalBusIndex = o_idx;
    return true;
  }
}

/**
 * Set the global index of the bus
 * @param idx: local index of bus
 * @param g_idx: global index to be assigned to bus
 * @return: false if no bus exists for idx
 */
bool setGlobalBusIndex(int idx, int g_idx)
{
  if (idx < 0 || idx >= p_buses.size()) {
    return false;
  } else {
    p_buses[idx]->p_globalBusIndex = g_idx;
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
  if (idx < 0 || idx >= p_branches.size()) {
    return false;
  } else {
    p_branches[idx]->p_globalBranchIndex = g_idx;
    return true;
  }
}

/**
 * Set the original index of the bus at the "from" end of branch
 * @param idx: local index of branch
 * @param b_idx: original index of "from" bus for this branch
 * @return: false if no branch exists for idx
 */
bool setOriginalBusIndex1(int idx, int b_idx)
{
  if (idx < 0 || idx >= p_branches.size()) {
    return false;
  } else {
    p_branches[idx]->p_originalBusIndex1 = b_idx;
    return true;
  }
}

/**
 * Set the original index of the bus at the "to" end of branch
 * @param idx: local index of branch
 * @param b_idx: original index of "to" bus for this branch
 * @return: false if no branch exists for idx
 */
bool setOriginalBusIndex2(int idx, int b_idx)
{
  if (idx < 0 || idx >= p_branches.size()) {
    return false;
  } else {
    p_branches[idx]->p_originalBusIndex2 = b_idx;
    return true;
  }
}

/**
 * Set the global index of the bus at the "from" end of branch
 * @param idx: local index of branch
 * @param b_idx: global index of "from" bus for this branch
 * @return: false if no branch exists for idx
 */
bool setGlobalBusIndex1(int idx, int b_idx)
{
  if (idx < 0 || idx >= p_branches.size()) {
    return false;
  } else {
    p_branches[idx]->p_globalBusIndex1 = b_idx;
    return true;
  }
}

/**
 * Set the global index of the bus at the "to" end of branch
 * @param idx: local index of branch
 * @param b_idx: global index of "to" bus for this branch
 * @return: false if no branch exists for idx
 */
bool setGlobalBusIndex2(int idx, int b_idx)
{
  if (idx < 0 || idx >= p_branches.size()) {
    return false;
  } else {
    p_branches[idx]->p_globalBusIndex2 = b_idx;
    return true;
  }
}

/**
 * Set the local index of the bus at the "from" end of branch
 * @param idx: local index of branch
 * @param b_idx: local index of "from" bus for this branch
 * @return: false if no branch exists for idx
 */
bool setLocalBusIndex1(int idx, int b_idx)
{
  if (idx < 0 || idx >= p_branches.size()) {
    return false;
  } else {
    p_branches[idx]->p_localBusIndex1 = b_idx;
    return true;
  }
}

/**
 * Set the local index of the bus at the "to" end of branch
 * @param idx: local index of branch
 * @param b_idx: local index of "to" bus for this branch
 * @return: false if no branch exists for idx
 */
bool setLocalBusIndex2(int idx, int b_idx)
{
  if (idx < 0 || idx >= p_branches.size()) {
    return false;
  } else {
    p_branches[idx]->p_localBusIndex2 = b_idx;
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
  if (idx < 0 || idx >= p_buses.size()) {
    return false;
  } else {
    p_buses[idx]->p_activeBus = flag;
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
  if (idx < 0 || idx >= p_branches.size()) {
    return false;
  } else {
    p_branches[idx]->p_activeBranch = flag;
    return true;
  }
}

/**
 * Clear the list of neighbors for the bus at idx
 * @param idx: local index of bus
 * @return: false if no bus exists for idx
 */
bool clearBranchNeighbors(int idx)
{
  if (idx < 0 || idx >= p_buses.size()) {
    return false;
  } else {
    p_buses[idx]->p_branchNeighbors.clear();
    return true;
  }
}

/**
 * Add local index for a branch attached to bus at idx
 * @param idx: local index of bus
 * @return: false if no bus exists for idx
 */
bool addBranchNeighbor(int idx, int br_idx)
{
  if (idx < 0 || idx >= p_buses.size()) {
    return false;
  } else {
    p_buses[idx]->p_branchNeighbors.push_back(br_idx);
    return true;
  }
}

// Bus and Branch accessors

/**
 * Get status of the bus (local or ghosted)
 * @param idx: local index of bus
 * @return: true if bus is locally held, false if it is ghosted 
 */
bool getActiveBus(int idx)
{
  if (idx >= 0 && idx < p_buses.size()) {
    return p_buses[idx]->p_activeBus;
  } else {
    printf("gridpack::network::getActiveBus: illegal index: %d size: %d\n",
            idx,p_buses.size());
    // TODO: some kind of error
  }
}

/**
 * Get original index of the bus
 * @param idx: local index of bus
 * @return: original index of bus 
 */
int getOriginalBusIndex(int idx)
{
  if (idx >= 0 && idx < p_buses.size()) {
    return p_buses[idx]->p_originalBusIndex;
  } else {
    // TODO: some kind of error
  }
}

/**
 * Get global index of the bus
 * @param idx: local index of bus
 * @return: global index of bus 
 */
int getGlobalBusIndex(int idx)
{
  if (idx >= 0 && idx < p_buses.size()) {
    return p_buses[idx]->p_globalBusIndex;
  } else {
    // TODO: some kind of error
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
    return p_buses[idx]->p_bus;
  }
}

/**
 * Get status of the branch (local or ghosted)
 * @param idx: local index of branch
 * @return: true if branch is locally held, false if it is ghosted 
 */
bool getActiveBranch(int idx)
{
  if (idx >= 0 && idx < p_branches.size()) {
    return p_branches[idx]->p_activeBranch;
  } else {
    // TODO: some kind of error
  }
}

/**
 * Retrieve a pointer to an existing branch
 * @param idx: local index of requested branch
 * @return: a pointer to the requested branch
 */
boost::shared_ptr<_branch> getBranch(int idx)
{
  if (idx<0 || idx >= p_branches.size()) {
    // TODO: Some kind of error
  } else {
    return p_branches[idx].p_branch;
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
    return p_buses[idx]->p_data;
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
    return p_branches[idx].p_data;
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
  return p_buses[idx]->p_branchNeighbors;
}

/**
 * Return list of buses connected to central bus via one branch
 * @param bus: local bus index
 * @return: vector of local bus indices
 */
std::vector<int> getConnectedBuses(int idx) const
{
  std::vector<int> branches = p_buses[idx]->p_branchNeighbors;
  int size = branches.size();
  std::vector<int> ret;
  int i, j;
  for (i=0; i<size; i++) {
    j = branches[i];
    if (p_branches[j]->p_localBusIndex1 != idx) {
      ret.push_back(p_branches[j]->p_localBusIndex1);
    } else {
      ret.push_back(p_branches[j]->p_localBusIndex2);
    }
  }
  return ret;
}

/**
 * Return indices of buses at either end of branch
 * @param bus: local branch index
 * @param bus1: local index of bus at one end of branch
 * @param bus2: local index of bus at other end of branch
 * @return: false if branch not found
 */
void getBranchEndpoints(int idx, int *bus1, int *bus2) const
{
  if (idx<0 || idx >= p_branches.size()) {
    // TODO: some kind of error
  } else {
    *bus1 = p_branches[idx]->p_localBusIndex1;
    *bus2 = p_branches[idx]->p_localBusIndex2;
  }
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
  int size = p_branches.size();
  int new_id = 0;
  for (i=0; i<size; i++) {
    if (p_branches[i]->p_activeBranch) {
      p_branches[new_id] = p_branches[i];
      branches.insert(std::pair<int, int>(i,new_id));
      new_id++;
    }
  }
  // clean up the ends of the branch vectors
  size = size - new_id;
  for (i=0; i<size; i++) {
    p_branches.pop_back();
  }

  // remove inactive buses
  size = p_buses.size();
  new_id = 0;
  for (i=0; i<size; i++) {
    if (p_buses[i]->p_activeBus) {
      p_buses[new_id] = p_buses[i];
      buses.insert(std::pair<int, int>(i,new_id));
      new_id++;
    }
  }
  // clean up the ends of the bus vectors
  size = size - new_id;
  for (i=0; i<size; i++) {
    p_buses.pop_back();
  }

  // reset all local indices
  size = p_branches.size();
  for (i=0; i<size; i++) {
    p = buses.find(p_branches[i]->p_localBusIndex1);
    if (p != buses.end()) {
      p_branches[i]->p_localBusIndex1 = p->second;
    } else {
      p_branches[i]->p_localBusIndex1 = -1;
    }
    p = buses.find(p_branches[i]->p_localBusIndex2);
    if (p != buses.end()) {
      p_branches[i]->p_localBusIndex2 = p->second;
    } else {
      p_branches[i]->p_localBusIndex2 = -1;
    }
  }

  int jsize;
  size = p_buses.size();
  for (i=0; i<size; i++) {
    std::vector<int> neighbors = p_buses[i]->p_branchNeighbors;
    jsize = neighbors.size();
    p_buses[i]->p_branchNeighbors.clear();
    for (j=0; j<jsize; j++) {
      p = branches.find(neighbors[j]);
      if (p != branches.end()) {
        p_buses[i]->p_branchNeighbors.push_back(p->second);
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
   * Vector of bus data and objects
   */
  std::vector<boost::shared_ptr<gridpack::network::BusData<_bus> > > p_buses;

  /**
   * Vector of branch data and objects
   */
  std::vector<boost::shared_ptr<gridpack::network::BranchData<_branch> > > p_branches;

  /**
   * Parallel environment for network
   */
#if 0
  ParallelEnv p_configuration;
#endif

  int p_refBus;
};
}  //namespace network
}  //namespace gridpack

#endif
