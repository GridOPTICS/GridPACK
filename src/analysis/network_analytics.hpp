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
 * @date   2024-05-15 09:12:00 d3g096
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
#include <tuple>
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
   * Get the total number of generators on a specific bus in the network
   * @return total number of generators on the bus
   */
  int numGenerators(const int& bus_idx)
  {
    int result(0);

    int lbusidx = p_network->getOriginalBusIndex(bus_idx);

    if (lbusidx >= 0) {
      if (p_network->getActiveBus(lbusidx)) {
        result += p_network->getBus(lbusidx)->numGenerators();
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
   * Get the total number of loads on a specific bus in the network
   * @return total number of loads on the bus
   */
  int numLoads(const int& bus_idx)
  {
    int result(0);

    int lbusidx = p_network->getOriginalBusIndex(bus_idx);

    if (lbusidx >= 0) {
      if (p_network->getActiveBus(lbusidx)) {
        result += p_network->getBus(lbusidx)->numLoads();
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
   * Get the total number of storage units on a specific bus in the network
   * @return total number of storage units on the bus
   */
  int numStorage(const int& bus_idx)
  {
    int result(0);

    int lbusidx = p_network->getOriginalBusIndex(bus_idx);

    if (lbusidx >= 0) {
      if (p_network->getActiveBus(lbusidx)) {
        result += p_network->getBus(lbusidx)->numStorage();
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

  /**
   * Get the total number of lines in a specific branch of the network
   * @return total number of lines in branch
   */
  int numLines(const int& branch_idx)
  {
    int result(0);
    int lbranchidx = p_network->getOriginalBranchIndex(branch_idx);

    if (lbranchidx >= 0) {
      if (p_network->getActiveBranch(lbranchidx)) {
        result += p_network->getBranch(lbranchidx)->numLines();
      }
    }
    p_network->communicator().sum(&result, 1);

    return result;
  };


  /// Get the number of buses
  /**
   * @e not collective
   *
   * @return number of buses in entire network
   **/
  int totalBuses(void)
  {
    return p_network->totalBuses();
  }

  /// Get the number of branches
  /**
   * @e not collective
   *
   * @return number of branches in entire network
   **/
  int totalBranches(void)
  {
    return p_network->totalBranches();
  }

  /// Get the branch indexes connected to a bus 
  /**
   * @e collective
   *
   * @param idx original bus index
   *
   * @result vector of original branch indexes
   **/
  std::vector<int> getConnectedBranches(const int& oidx) const
  {
    int rsize(0);
    std::vector<int> result;
    std::vector<int> lbusi(p_network->getLocalBusIndices(oidx));

    // look for the active bus index -- it should only be active on
    // one process -- result should only be filled on one process

    if (!lbusi.empty()) {
      for (auto ib : lbusi) {
        if (p_network->getActiveBus(ib)) {
          typename NetworkType::BusPtr bus(p_network->getBus(ib));
          std::vector<int> lbranchi(p_network->getConnectedBranches(ib));
          for (auto ibr : lbranchi) {
            result.push_back(p_network->getGlobalBranchIndex(ibr));
          }
          rsize = result.size();
        }
      }
    }

    p_network->communicator().sum(&rsize, 1);

    // only one process should have non-zeros in result

    if (result.empty()) {
      result.resize(rsize, 0);
    }
    
    p_network->communicator().sum(&result[0], rsize);

    return result;
  }

    
  /// Get the bus indexes connected to a branch
  /**
   * @e collective
   *
   * @param idx original branch index
   * @param fbus pointer to location to put from-bus (original) index
   * @param tbus pointer to location to put to-bus (original) index
   **/
  void getBranchEndpoints(const int& idx, int *fbus, int *tbus) const
  {
    int fresult(0), tresult(0);
    int lbranchidx = p_network->getOriginalBranchIndex(idx);

    if (lbranchidx >= 0) {
      if (p_network->getActiveBranch(lbranchidx)) {
        p_network->getOriginalBranchEndpoints(lbranchidx, &fresult, &tresult);
      }
    }
    p_network->communicator().sum(&fresult, 1);
    p_network->communicator().sum(&tresult, 1);

    *fbus = fresult;
    *tbus = tresult;
  }
  

  private:

  NetworkPtr p_network;
};


} // namespace gridpack
} // namespace analysis

#endif
