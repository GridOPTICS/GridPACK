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
   * Clean all ghost buses and branches from the system. This can be used
   * before repartitioning the network
   */
   void Clean(void);

   /**
    * Update the ghost values of this field. This is a
    * collective operation across all processors.
    * @param field: pointer to the network field that must be
    *       updated
    */
   void UpdateField(BaseField *field);

   // TODO: Need operations to move buses and branches to other processors in partitioners

protected:

   /**
    * Protected copy constructor to avoid unwanted copies.
    */
   BaseNetwork(const BaseNetwork& old);

private:

   vector<bool> p_activeBus;

   vector<bool> p_activeBranch;

   BusField<int> p_originalIndex;

   BusField<int> p_globalIndex;

   BusField<std::vector<int> > p_busNeighbors;

   BusField<std::vector<int> > p_branchNeighbors;

   BranchField<int> p_globalBranchIndex1;

   BranchField<int> p_globalBranchIndex2;

   BranchField<int> p_localBranchIndex1;

   BranchField<int> p_localBranchIndex2;

   std::map<std::string, BusField*> p_busFields;

   std::map<std::string, BranchField*> p_branchFields;

};

#endif
