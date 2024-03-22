// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   network_analytics.hpp
 * @author Bruce Palmer
 * @date   2024-03-22
 * 
 * @brief  
 * This is a utility that can be used to extract properties of the network
 * using the network analytics interface for the individual network components.
 * These include properties such as the total number of generators, loads,
 * and lines in a system as well as any other properties that may be of
 * interest.
 * 
 */

// -------------------------------------------------------------

#ifndef _network_analytics_hpp_
#define _network_analytics_hpp_

#include <ga.h>
#include <map>
#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parallel/communicator.hpp"
#include "gridpack/network/base_network.hpp"

namespace gridpack {
namespace analysis {

template <class _network> class NetworkAnalytics {
  public:
  typedef _network NetworkType;
  typedef boost::shared_ptr<NetworkType> NetworkPtr;

  /**
   * Constructor
   * @param network network that is used for reporting analytics
   */
  NetworkAnalytics(NetworkPtr network): p_network(network)
  {
    /* initialize network analytics functionality */
    int nbus = p_network->numBuses();
    for (size_t i=0; i<nbus; i++) {
      p_network->getBus(i)->setData(p_network->getBusData(i));
    }
    int nbranch = p_network->numBranches();
    for (size_t i=0; i<nbranch; i++) {
      p_network->getBranch(i)->setData(p_network->getBranchData(i));
    }
  }

  /**
   * Destructor
   */
  virtual ~NetworkAnalytics(void)
  {
  }

  /**
   * Get the total number of generators in the network
   * @return total number of generators
   */
  int numGenerators()
  {
    int result(0);
    const int nBus(p_network->numBuses());

    for (int i = 0; i < nBus; ++i) {
      if (p_network->getActiveBus(i)) {
        result += p_network->getBus(i)->numGenerators();
      }
    }
    p_network->communicator().sum(&result, 1);
    return result;
  }

  /**
   * Get the total number of loads in the network
   * @return total number of loads
   */
  int numLoads()
  {
    int result(0);
    const int nBus(p_network->numBuses());

    for (int i = 0; i < nBus; ++i) {
      if (p_network->getActiveBus(i)) {
        result += p_network->getBus(i)->numLoads();
      }
    }
    p_network->communicator().sum(&result, 1);
    return result;
  }

  /**
   * Get the total number of storage units in the network
   * @return total number of storage units
   */
  int numStorage()
  {
    int result(0);
    const int nBus(p_network->numBuses());

    for (int i = 0; i < nBus; ++i) {
      if (p_network->getActiveBus(i)) {
        result += p_network->getBus(i)->numStorage();
      }
    }
    p_network->communicator().sum(&result, 1);
    return result;
  }

  /**
   * Get the total number of lines in the network
   * @return total number of lines
   */
  int numLines()
  {
    int result(0);
    const int nBranch(p_network->numBranches());

    for (int i = 0; i < nBranch; ++i) {
      if (p_network->getActiveBranch(i)) {
        result += p_network->getBranch(i)->numLines();
      }
    }
    p_network->communicator().sum(&result, 1);
    return result;
  };

  private:

  NetworkPtr p_network;
};


} // namespace gridpack
} // namespace analysis

#endif
