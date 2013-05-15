// -------------------------------------------------------------
/**
 * @file   base_network.cpp
 * @author Bruce Palmer, William Perkins
 * @date   April 3, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <map>
#include "gridpack/parallel/distribution.hpp"
#include "gridpack/network/base_network.hpp"

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
/**
 * Default constructor.
 */
BaseNetwork::BaseNetwork(void)
{
  //TODO: Get default parallel configuration for world group
  p_refBus = -1;
}

/**
 * Constructor with parallel configuration for running on
 * subgroups or communicators other than the world
 * communicator
 */
BaseNetwork::BaseNetwork(ParallelEnv configuration)
{
  p_configuration = configuration;
  p_refBus = -1;
}

/**
 * Default destructor.
 */
BaseNetwork::~BaseNetwork(void)
{
}

/**
 * Add a bus locally to the network
 * @param idx: original index of  bus
 */
void BaseNetwork::addBus(int idx)
{
  p_originalBusIndex.push_back(idx);
  p_activeBus.push_back(true);
  std::map<std::string, stlplus::smart_ptr<BusField> >::iterator bus;
  bus = p_busFields.begin();
  while (bus != p_busFields.end()) {
    bus->second->append();
    bus++;
  }
}

/**
 * Add a branch locally to the network. A branch is defined by
 * buses at either end
 * @param idx1: original bus index of bus 1
 * @param idx2: original bus index of bus 2
 */
void BaseNetwork::addBranch(int idx1, int idx2)
{
  p_originalBranchIndex1.push_back(idx1);
  p_originalBranchIndex2.push_back(idx2);
  p_activeBranch.push_back(true);
  std::map<std::string, stlplus::smart_ptr<BranchField> >::iterator branch;
  branch = p_branchFields.begin();
  while (branch != p_branchFields.end()) {
    branch->second->append;
    branch++;
  }
}

/**
 * Designate a bus as a reference bus.
 * @param idx: local index of bus
 */
void BaseNetwork::setReferenceBus(int idx)
{
  p_refBus = idx;
}

/**
 * Return index of reference bus.
 * @return: local index of reference bus. If reference bus is not on this
 * processor then return -1.
 */
int BaseNetwork::getReferenceBus(void) const
{
  return p_refBus;
}

/**
 * Add a new field to the network buses
 * @param name: a string designating the name of the field
 * @param field: a pointer to the BusField being added to the
 *       network
 */
void BaseNetwork::addBusField(std::string name, stlplus::smart_ptr<BusField> field)
{
  // check size of new field to see if it is too large
  if (field->size() > p_activeBus.size()) {
    //TODO: Some kind of failure
  } else if (field->size()<p_activeBus.size()) {
    // pad out field so it matches size of existing fields
    int oldsize = p_activeBus.size();
    int newsize = field->size();
    int i;
    for (i=0; i<oldsize-newsize; i++) {
      field->append();
    }
  }
  p_busFields.insert(std::pair<std::string, stlplus::smart_ptr<BusField> >(name,field));
}

/**
 * Add a new field to the network branches
 * @param name: a string designating the name of the field
 * @param field: a pointer to the BranchField being added to
 *       the network
 */
void BaseNetwork::addBranchField(std::string name, stlplus::smart_ptr<BranchField> field)
{
  // check size of new field to see if it is too large
  if (field->size() > p_activeBus.size()) {
    //TODO: Some kind of failure
  } else if (field->size()<p_activeBranch.size()) {
    // pad out field so it matches size of existing fields
    int oldsize = p_activeBranch.size();
    int newsize = field->size();
    int i;
    for (i=0; i<oldsize-newsize; i++) {
      field->append();
    }
  }
  p_branchFields.insert(std::pair<std::string, stlplus::smart_ptr<BranchField> >(name,field));
}

/**
 * Retrieve a pointer to an existing bus field
 * @param name: a string representing the name of the desired
 *       field
 * @return: a pointer to the requested field. If the field is
 *       not found, the pointer is null
 */
smart_ptr<BusField> BaseNetwork::getBusField(std::string name) const
{
  std::map<std::string, stlplus::smart_ptr<BusField> >::iterator bus;
  bus = p_busFields.find(name);
  if (bus != p_busFields.end()) {
    return bus->second;
  } else {
    return (BusField*)0;
  }
}

/**
 * Retrieve a pointer to an existing branch field
 * @param name: a string representing the name of the desired
 *       field
 * @return: a pointer to the requested field. If the field is
 *       not found, the pointer is null
 */
smart_ptr<BranchField> BaseNetwork::getBranchField(std::string name) const
{
  std::map<std::string, stlplus::smart_ptr<BranchField> >::iterator branch;
  branch = p_branchFields.find(name);
  if (branch != p_branchFields.end()) {
    return branch->second;
  } else {
    return (BranchField*)0;
  }
}

/**
 * Delete an existing bus field
 * @param name: a string representing the name of the field
 *       to be deleted
 */
void BaseNetwork::deleteBusField(std::string name)
{
  std::map<std::string, stlplus::smart_ptr<BusField> >::iterator bus;
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
void BaseNetwork::deleteBranchField(std::string name)
{
  std::map<std::string, stlplus::smart_ptr<BranchField> >::iterator branch;
  branch = p_busFields.find(name);
  if (branch != p_branchFields.end()) {
    p_branchFields.erase(branch);
  }
}

/**
 * Return list of branches connected to bus
 * @param bus: local bus index
 * @return: vector of local branch indices
 */
vector<int> BaseNetwork::getConnectedBranches(int idx) const
{
  return p_branchNeighbors[idx];
}

/**
 * Return list of buses connected to central bus via one branch
 * @param bus: local bus index
 * @return: vector of local bus indices
 */
vector<int> BaseNetwork::getConnectedBuses(int idx) const
{
  vector<int> branches = p_branchNeighbors[idx];
  int size = branches.size();
  vector<int> ret;
  int i, j;
  for (i=0; i<size; i++) {
    j = branches[i];
    if (p_localBranchIndex1[j] != p_idx) {
      ret.push_back(p_localBranchIndex1[j]);
    } else {
      ret.push_back(p_localBranchIndex2[j]);
    }
  }
}

/**
 * Return indices of buses at either end of branch
 * @param bus: local branch index
 * @param bus1: local index of bus at one end of branch
 * @param bus2: local index of bus at other end of branch
 */
void BaseNetwork::getBusEndpoints(int idx, int *bus1, int *bus2) const
{
  *bus1 = p_localBranchIndex1[idx];
  *bus2 = p_localBranchIndex2[idx];
}

/**
 * Clean all ghost buses and branches from the system. This can be used
 * before repartitioning the network
 */
void BaseNetwork::clean(void)
{
  std::map<int new_id, int old_id> buses;
  std::map<int new_id, int old_id> branches;
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
      branches.insert(std::pair<int, int>(i,new_id));
      std::map<std::string, stlplus::smart_ptr<BranchField> >::iterator branch;
      branch = p_branchFields->begin();
      while(branch !=  p_branchFields->end()) {
        (*(branch->second))[new_id] = (*(branch->second))[i];
      }
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
    std::map<std::string, stlplus::smart_ptr<BranchField> >::iterator branch;
    branch = p_branchFields->begin();
    while(branch !=  p_branchFields->end()) {
      branch->second->pop_back();
    }
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
      buses.insert(std::pair<int, int>(i,new_id));
      std::map<std::string, stlplus::smart_ptr<BranchField> >::iterator bus;
      bus = p_busFields->begin();
      while(bus !=  p_busFields->end()) {
        (*(bus->second))[new_id] = (*(bus->second))[i];
      }
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
    std::map<std::string, stlplus::smart_ptr<BranchField> >::iterator bus;
    bus = p_busFields->begin();
    while(bus !=  p_busFields->end()) {
      bus->second->pop_back();
    }
  }
  // reset all local indices
  size = p_activeBranch.size();
  for (i=0; i<size; i++) {
    p_localBranchIndex1[i] = buses[p_localBranchIndex1[i]];
    p_localBranchIndex2[i] = buses[p_localBranchIndex2[i]];
    std::map<std::string, stlplus::smart_ptr<BranchField> >::iterator branch;
    branch = p_branchFields->begin();
    while(branch !=  p_branchFields->end()) {
      branch->second->setIndex(branches[branch->second->getIndex()]);
    }
  }
  int jsize;
  size = p_activeBus.size();
  for (i=0; i<size; i++) {
    jsize = p_branchNeighbors[i].size();
    for (j=0; j<jsize; j++) {
      (p_branchNeighbors[i])[j] = branches[(p_branchNeighbors[i])[j]];
    }
    jsize = p_busNeighbors[i].size();
    for (j=0; j<jsize; j++) {
      (p_busNeighbors[i])[j] = buses[(p_busNeighbors[i])[j]];
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
void BaseNetwork::updateField(char *field)
{
}

/**
 * Protected copy constructor to avoid unwanted copies.
 */
BaseNetwork::BaseNetwork(const BaseNetwork& old)
{
}
