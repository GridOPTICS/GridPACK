// -------------------------------------------------------------
/**
 * @file   fields.cpp
 * @author Bruce Palmer
 * @date   April 25, 2013
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
#include "gridpack/network/fields.hpp"
// -------------------------------------------------------------
//  class BaseField:
//  This class implements some basic functions that can be
//  expected from any field on the network.
// -------------------------------------------------------------

/**
 * Constructor
 * @param network: Network associated with field
 */
BaseField::BaseField(BaseNetwork *network)
{
  p_network = network;
  p_type = elem;
}

/**
 * Destructor
 */
BaseField::~BaseField(void)
{
  p_vector.clear();
}

/**
 * Index operator
 * @param index: index of element in field
 * @return: field element at local index
 */
elem& BaseField::operator[] (int index)
{
  return p_vector[i];
}

/**
 * Return size of field
 */
int BaseField::size(void)
{
  return p_vector.size();
}

/**
 * Return whether the grid element at the index location is
 * active (local) or inactive (ghost)
 * @param index: index of element in field
 */
bool BaseField::active(int index)
{
  return p_network->p_activeBus[index];
}

/**
 * Add another element to the field
 * @param new_elem: new element to be added to field
 * @return: index of new element
 */
int BaseField::append(elem new_elem)
{
  p_vector.push_back(elem); 
  return p_vector.size()-1;
}

/**
 * Add another element to the field. Use default element
 * @return: index of new element
 */
int BaseField::append(void)
{
  p_vector.push_back(p_type); 
  return p_vector.size()-1;
}

/**
 * Delete element from field
 * @param index: index of element to be deleted
 * @return: success or failure of delete operation
 */
bool BaseField::delete(int index)
{
  int size = p_vector.size();
  if (index >= size || index < 0) {
    return false;
  }
  if (index == size-1) {
    p_vector.pop_back();
    return true;
  }
  vector<elem>::iterator p;
  p = vector.begin();
  for (i=0; i<index; i++) p++;
  p_vector.erase(p);
  return true;
}

/**
 * Clear all elements from field
 */
void BaseField::clear(void)
{
  p_vector.clear();
}

// -------------------------------------------------------------
//  class BusField:
//  This class represents fields defined on the network buses
// -------------------------------------------------------------
/**
 * Constructor
 * @param network: Network associated with field
 */
BusField::BusField(BaseNetwork *network)
{
  p_network = network;
  p_type = elem;
}

/**
 * Destructor
 */
BusField::~BusField(void)
{
  p_vector.clear();
}

/**
 * Get the local indices of the branches that are attached
 * to the bus element indicated by index
 * @param index: index of bus element
 * @return: list of local indices of branches attached to
 *         this bus
 */
vector<int> BusField::neighborBranches(int index);
{
  return network->p_branchNeighbors[index];
}

/**
 * Get the local indices of the buses that are connected
 * to this bus via a branch element
 * @param index: index of bus element
 * @return: list of local indices of neighboring buses
 */
vector<int> BusField::neighborBuses(int index);
{
  return network->p_busNeighbors[index];
}

// -------------------------------------------------------------
//  class BusField:
//  This class represents fields defined on the network branches
// -------------------------------------------------------------
/**
 * Constructor
 * @param network: Network associated with field
 */
BranchField::BranchField(BaseNetwork *network)
{
  p_network = network;
  p_type = elem;
}

/**
 * Destructor
 */
BranchField::~BranchField(void)
{
  p_vector.clear();
}

/**
 * Get one of the terminal buses on this branch
 * @param index: index of branch element
 * @return: index to the bus at one end of the branch
 */
int BranchField::getBus1(int index)
{
  return network->p_localBranchIndex1[index];
}

/**
 * Get the other terminal bus on this branch
 * @param index: index of branch element
 * @return: index to the bus at the other end of the
 *        branch
 */
int BranchField::getBus2(int index)
{
  return network->p_localBranchIndex2[index];
}
