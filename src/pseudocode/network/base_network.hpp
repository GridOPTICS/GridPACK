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
#include "gridpack/parallel/distribution.hpp"
// -------------------------------------------------------------
//  class BaseField:
//  This class implements some basic functions that can be
//  expected from any field on the network.
// -------------------------------------------------------------
template <class elem>
class BaseField  {
  public:
    /**
     * Constructor
     */
    BaseField(void);

    /**
     * Destructor
     */
    ~BaseField(void);

    /**
     * Index operator
     * @param index: index of element in field
     */
    elem& operator[] (int index);

    /**
     * Return size of field
     */
    int Size(void);

    /**
     * Return whether the grid element at the index location is
     * active (local) or inactive (ghost)
     * @param index: index of element in field
     */
    bool Active(int index);

    friend class BaseNetwork;
  private:
    /**
     * Add another element to the field
     * @return: index of new element
     */
    int Append(elem *new_elem);

    /**
     * Delete element from field
     * @param index: index of element to be deleted
     * @return: success or failure of delete operation
     */
    bool Delete(int index);

    /**
     * Clear all elements from field
     */
    void Clear(void);

};

// -------------------------------------------------------------
//  class BusField:
//  This class represents fields defined on the network buses
// -------------------------------------------------------------
template <class elem>
class BusField : public BaseField<elem> {
  public:
    /**
     * Constructor
     */
    BusField(void);

    /**
     * Destructor
     */
    ~BusField(void);

    /**
     * Get the local indices of the branches that are attached
     * to the bus element indicated by index
     * @param index: index of bus element
     * @return: list of local indices of branches attached to
     *         this bus
     */
    vector<int> NeighborBranches(int index);

    /**
     * Get the local indices of the buses that are connected
     * to this bus via a branch element
     * @param index: index of bus element
     * @return: list of local indices of neighboring buses
     */
    vector<int> NeighborBuses(int index);
};

// -------------------------------------------------------------
//  class BusField:
//  This class represents fields defined on the network branches
// -------------------------------------------------------------
template <class elem>
class BranchField : public BaseField<elem> {
  public:
    /**
     * Constructor
     */
    BaseField(void);

    /**
     * Destructor
     */
    ~BaseField(void);

    /**
     * Get one of the terminal buses on this branch
     * @param index: index of branch element
     * @return: index to the bus at one end of the branch
     */
    int GetBus1(int index);

    /**
     * Get the other terminal bus on this branch
     * @param index: index of branch element
     * @return: index to the bus at the other end of the
     *        branch
     */
    int GetBus2(int index);

};

// -------------------------------------------------------------
//  class BaseNetwork:
//  This is the base class for creating distributed networks. It
//  is basically something that supports the network topology,
//  can have fields added and subtracted to the buses (nodes)
//  and branches (edges) of the network, and implements ghost
//  updates on the network. The BaseNetwork class does not
//  contain the partitioner but it does contain methods that
//  allow the partitioner to move grid elements around and
//  create ghost elements.
// -------------------------------------------------------------
class BaseNetwork {
public:
  /**
   * The field classes will need to be able to extract information
   * about the topology of the network directly from the network
   */
  friend class BaseField;
  friend class BusField;
  friend class BranchField;

  /**
   * Default constructor.
   */
  BaseNetwork(void);

  /**
   * Constructor with parallel configuration for running on
   * subgroups or communicators other than the world
   * communicator
   */
  BaseNetwork(ParallelEnv configuration);

  /**
   * Default destructor.
   */
  ~BaseNetwork(void);

  /**
   * Add a bus locally to the network
   * @param idx: global index of  bus
   */
  void AddBus(int idx);

  /**
   * Add a branch locally to the network. A branch is defined by
   * buses at either end
   * @param idx1: global bus index of bus 1
   * @param idx2: global bus index of bus 2
   */
  void AddBranch(int idx1, int idx2);

  /**
   * Designate a bus as a reference bus.
   * @param idx: global index of bus
   */
  void SetReferenceBus(int idx);

  /**
   * Add a new field to the network buses
   * @param name: a string designating the name of the field
   * @param field: a pointer to the BusField being added to the
   *       network
   */
  void AddBusField(std::string name, BusField *field)

  /**
   * Add a new field to the network branches
   * @param name: a string designating the name of the field
   * @param field: a pointer to the BranchField being added to
   *       the network
   */
  void AddBranchField(std::string name, BranchField *field)

  /**
   * Retrieve a pointer to an existing bus field
   * @param name: a string representing the name of the desired
   *       field
   * @return: a pointer to the requested field. If the field is
   *       not found, the pointer is null
   */
   BusField* GetBusField(std::string name);

  /**
   * Retrieve a pointer to an existing branch field
   * @param name: a string representing the name of the desired
   *       field
   * @return: a pointer to the requested field. If the field is
   *       not found, the pointer is null
   */
   BranchField* GetBranchField(std::string name);

  /**
   * Delete an existing bus field
   * @param name: a string representing the name of the field
   *       to be deleted
   */
   void DeleteBusField(std::string name);

  /**
   * Delete an existing branch field
   * @param name: a string representing the name of the field
   *       to be deleted
   */
   void DeleteBranchField(std::string name);

   /**
    * Update the ghost values of this field. This is a
    * collective operation across all processors.
    * @param field: pointer to the network field that must be
    *       updated
    */
   void UpdateField(BaseField *field);

   // Need to move buses and branches to other processors in partitioners

protected:

   /**
    * Protected copy constructor to avoid unwanted copies.
    */
   BaseNetwork(const BaseNetwork& old);

private:

   BusField<bool> p_activeBus;

   BranchField<bool> p_activeBranch;

   BusField<int> p_originalIndex;

   BusField<int> p_globalIndex;

   BusField<std::vector<int> > p_busNeighbors;

   BusField<std::vector<int> > p_branchNeighbors;

   BranchField<int> p_globalBranchIndex1;

   BranchField<int> p_globalBranchIndex2;

   BranchField<int> p_localBranchIndex1;

   BranchField<int> p_localBranchIndex2;

   std::map<std:string, BusField*> p_busFields;

   std::map<std:string, BranchField*> p_branchFields;

};

#endif
