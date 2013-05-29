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
#include "gridpack/network/fields.hpp"
#include "gridpack/base_component/base_network_component.hpp"

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
class BaseNetwork {
public:
  /**
   * The field classes will need to be able to extract information
   * about the topology of the network directly from the network
   */
#if 0
  friend class BaseField;
  friend class BusField;
  friend class BranchField;
#endif

  /**
   * Default constructor.
   */
  BaseNetwork(void);

  /**
   * Constructor with parallel configuration for running on
   * subgroups or communicators other than the world
   * communicator
   */
#if 0
  BaseNetwork(ParallelEnv configuration);
#endif // Hold off on this until we work out details of parallel environment

  /**
   * Default destructor.
   */
  ~BaseNetwork(void);

  /**
   * Add a bus locally to the network
   * @param idx: original index of  bus
   */
  void addBus(int idx);

  /**
   * Add a branch locally to the network. A branch is defined by
   * buses at either end
   * @param idx1: original bus index of bus 1
   * @param idx2: original bus index of bus 2
   */
  void addBranch(int idx1, int idx2);

  /**
   * Designate a bus as a reference bus.
   * @param idx: local index of bus
   */
  void setReferenceBus(int idx);

  /**
   * Return index of reference bus.
   * @return: local index of reference bus. If reference bus is not on this
   * processor then return -1.
   */
  int getReferenceBus(void) const;

  /**
   * Add a new field to the network buses
   * @param name: a string designating the name of the field
   * @param field: a pointer to the BusField being added to the
   *       network
   */
  void addBusField(std::string name,
       boost::shared_ptr<gridpack::network::BusField<gridpack::component::BaseNetworkComponent> > field);

  /**
   * Add a new field to the network branches
   * @param name: a string designating the name of the field
   * @param field: a pointer to the BranchField being added to
   *       the network
   */
  void addBranchField(std::string name,
       boost::shared_ptr<gridpack::network::BranchField<gridpack::component::BaseNetworkComponent> > field);

  /**
   * Retrieve a pointer to an existing bus field
   * @param name: a string representing the name of the desired
   *       field
   * @return: a pointer to the requested field. If the field is
   *       not found, the pointer is null
   */
  boost::shared_ptr<gridpack::network::BusField<gridpack::component::BaseNetworkComponent> > getBusField(std::string name);

  /**
   * Retrieve a pointer to an existing branch field
   * @param name: a string representing the name of the desired
   *       field
   * @return: a pointer to the requested field. If the field is
   *       not found, the pointer is null
   */
  boost::shared_ptr<gridpack::network::BranchField<gridpack::component::BaseNetworkComponent> > getBranchField(std::string name);

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
   * Return list of branches connected to bus
   * @param bus: local bus index
   * @return: vector of local branch indices
   */
  std::vector<int> getConnectedBranches(int idx) const;

  /**
   * Return list of buses connected to central bus via one branch
   * @param bus: local bus index
   * @return: vector of local bus indices
   */
  std::vector<int> getConnectedBuses(int idx) const;

  /**
   * Return indices of buses at either end of branch
   * @param bus: local branch index
   * @param bus1: local index of bus at one end of branch
   * @param bus2: local index of bus at other end of branch
   */
  void getBusEndpoints(int idx, int *bus1, int *bus2) const;

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
   * BusFields associated with buses. These can be accessed by name
   */
  std::map<std::string,
    boost::shared_ptr<gridpack::network::BusField<gridpack::component::BaseNetworkComponent> > > p_busFields;

  /**
   * BranchFields associated with buses. These can be accessed by name
   */
  std::map<std::string,
    boost::shared_ptr<gridpack::network::BranchField<gridpack::component::BaseNetworkComponent> > > p_branchFields;

  /**
   * Parallel environment for network
   */
#if 0
  ParallelEnv p_configuration;
#endif

  /**
   * Reference bus index
   */
  int p_refBus;
};
}  //namespace network
}  //namespace gridpack

#endif
