// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: network_topology_interface.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created February 13, 2024 by Perkins
// -------------------------------------------------------------


#ifndef _network_topology_hpp_
#define _network_topology_hpp_

#include <boost/shared_ptr.hpp>
#include <vector>
#include <gridpack/parallel/distributed.hpp>

namespace gridpack {
namespace network {

namespace gp = gridpack::parallel;

// -------------------------------------------------------------
// class NetworkTopologyInterface
//
// This is an abstract interface definition for querying the network
// topology.  This includes only bus/branch counts and connectivity
// -------------------------------------------------------------
class NetworkTopologyInterface
  : public gp::Distributed
{
public:

  /// Default constructor on world communicator
  NetworkTopologyInterface(void)
    : gp::Distributed(gp::Communicator())
  {};

  /// Constructor on a specific communicator
  explicit NetworkTopologyInterface(const gp::Communicator& comm)
    : gp::Distributed(comm)
  {};

  /// Destructor
  ~NetworkTopologyInterface(void)
  {};

  /// Get the number of buses
  /**
   * @e not collective
   *
   * @return number of buses in entire network
   **/
  virtual int totalBuses(void) = 0;

  /// Get the number of branches
  /**
   * @e not collective
   *
   * @return number of branches in entire network
   **/
  virtual int totalBranches(void) = 0;

  /// Get the branch indexes connected to a bus 
  /**
   * @e not collective
   *
   * @param idx local bus index
   *
   * @result vector of local branch indexes
   **/
  // virtual std::vector<int> getConnectedBranches(int idx) const = 0;

  /// Get the bus indexes connected to a branch
  /**
   * @e not collective
   *
   * @param idx local branch index
   * @param fbus pointer to location to put from-bus (local) index
   * @param tbus pointer to location to put to-bus (local) index
   **/
  // virtual void getBranchEndpoints(int idx, int *fbus, int *tbus) const = 0;

  

};

} // namespace gridpack
} // namespace network


#endif
