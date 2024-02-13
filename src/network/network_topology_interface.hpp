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
  virtual int totalBuses(void) = 0;

  /// Get the number of branches
  virtual int totalBranches(void) = 0;

};

} // namespace gridpack
} // namespace network


#endif
