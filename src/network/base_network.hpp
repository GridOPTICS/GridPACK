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
  void addBus(int idx);

  /**
   * Add a branch locally to the network. A branch is defined by
   * buses at either end
   * @param idx1: global bus index of bus 1
   * @param idx2: global bus index of bus 2
   */
  void addBranch(int idx1, int idx2);

  /**
   * Designate a bus as a reference bus.
   * @param idx: global index of bus
   */
  void setReferenceBus(int idx);

  /**
   * Add a new field to the network buses
   * @param name: a string designating the name of the field
   * @param field: a pointer to the BusField being added to the
   *       network
   */
  void addBusField(std::string name, BusField *field)

  /**
   * Add a new field to the network branches
   * @param name: a string designating the name of the field
   * @param field: a pointer to the BranchField being added to
   *       the network
   */
  void addBranchField(std::string name, BranchField *field)

  /**
   * Retrieve a pointer to an existing bus field
   * @param name: a string representing the name of the desired
   *       field
   * @return: a pointer to the requested field. If the field is
   *       not found, the pointer is null
   */
   BusField* getBusField(std::string name);

  /**
   * Retrieve a pointer to an existing branch field
   * @param name: a string representing the name of the desired
   *       field
   * @return: a pointer to the requested field. If the field is
   *       not found, the pointer is null
   */
   BranchField* getBranchField(std::string name);

  /**
   * Delete an existing bus field
   * @param name: a string representing the name of the field
   *       to be deleted
   */
   void deleteBusField(std::string name);

  /**
   * Delete an existing branch field
   * @param name: a string representing the name of the field
   *       to be deleted
   */
   void deleteBranchField(std::string name);

  /**
   * Clean all ghost buses and branches from the system. This can be used
   * before repartitioning the network
   */
   void clean(void);

   /**
    * Update the ghost values of this field. This is a
    * collective operation across all processors.
    * @param field: name of the network field that must be
    *       updated
    */
   void updateField(char *field);

   // TODO: Need operations to move buses and branches to other processors in partitioners

protected:

   /**
    * Protected copy constructor to avoid unwanted copies.
    */
   BaseNetwork(const BaseNetwork& old);

private:

   std::vector<bool> p_activeBus;

   std::vector<bool> p_activeBranch;

   std::vector<int> p_originalIndex;

   std::vector<int> p_globalIndex;

   std::vector<std::vector<int> > p_busNeighbors;

   std::vector<std::vector<int> > p_branchNeighbors;

   std::vector<int> p_globalBranchIndex1;

   std::vector<int> p_globalBranchIndex2;

   std::vector<int> p_localBranchIndex1;

   std::vector<int> p_localBranchIndex2;

   std::map<std::string, BusField*> p_busFields;

   std::map<std::string, BranchField*> p_branchFields;

};

#endif
