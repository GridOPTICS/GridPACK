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
  std::map<std::string, BusField*>::iterator bus;
  bus = p_busFields.begin();
  while (bus != p_busFields.end()) {
    delete bus->second;
    bus++;
  }
  std::map<std::string, BranchField*>::iterator branch;
  branch = p_branchFields.begin();
  while (branch != p_branchFields.end()) {
    delete branch->second;
    branch++;
  }
}

/**
 * Add a bus locally to the network
 * @param idx: global index of  bus
 */
void BaseNetwork::addBus(int idx)
{
  p_originalBusIndex.push_back(idx);
  p_activeBus.push_back(true);
  std::map<std::string, BusField*>::iterator bus;
  bus = p_busFields.begin();
  while (bus != p_busFields.end()) {
    bus->second->append();
    bus++;
  }
}

/**
 * Add a branch locally to the network. A branch is defined by
 * buses at either end
 * @param idx1: global bus index of bus 1
 * @param idx2: global bus index of bus 2
 */
void BaseNetwork::addBranch(int idx1, int idx2)
{
  p_originalBranchIndex1.push_back(idx1);
  p_originalBranchIndex2.push_back(idx2);
  p_activeBranch.push_back(true);
  std::map<std::string, BranchField*>::iterator branch;
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
int BaseNetwork::getReferenceBus(void)
{
  return p_refBus;
}

/**
 * Add a new field to the network buses
 * @param name: a string designating the name of the field
 * @param field: a pointer to the BusField being added to the
 *       network
 */
void BaseNetwork::addBusField(std::string name, BusField *field)
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
  p_busFields.insert(std::pair<std::string, *BusField>(name,field));
}

/**
 * Add a new field to the network branches
 * @param name: a string designating the name of the field
 * @param field: a pointer to the BranchField being added to
 *       the network
 */
void BaseNetwork::addBranchField(std::string name, BranchField *field)
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
  p_branchFields.insert(std::pair<std::string, *BranchField>(name,field));
}

/**
 * Retrieve a pointer to an existing bus field
 * @param name: a string representing the name of the desired
 *       field
 * @return: a pointer to the requested field. If the field is
 *       not found, the pointer is null
 */
BusField* BaseNetwork::getBusField(std::string name)
{
  std::map<std::string, BusField*>::iterator bus;
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
BranchField* BaseNetwork::getBranchField(std::string name)
{
  std::map<std::string, BranchField*>::iterator branch;
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
  std::map<std::string, BusField*>::iterator bus;
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
  std::map<std::string, BranchField*>::iterator branch;
  branch = p_busFields.find(name);
  if (branch != p_branchFields.end()) {
    p_branchFields.erase(branch);
  }
}

/**
 * Clean all ghost buses and branches from the system. This can be used
 * before repartitioning the network
 */
void BaseNetwork::clean(void)
{
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
