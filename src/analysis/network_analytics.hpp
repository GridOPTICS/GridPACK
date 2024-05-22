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

    int lidx(p_localFromGlobalBusIndex(bus_idx));

    if (lidx >= 0) {
      result = p_network->getBus(lidx)->numGenerators();
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
   * Get the number of loads on a specific bus in the network
   * @param bus_idx @e global bus index
   * @return total number of loads on the bus
   */
  int numLoads(const int& bus_idx)
  {
    int result(0);

    int lidx(p_localFromGlobalBusIndex(bus_idx));

    if (lidx >= 0) {
      result = p_network->getBus(lidx)->numLoads();
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
   * @param bus_idx @e global bus index
   * @return total number of storage units on the bus
   */
  int numStorage(const int& bus_idx)
  {
    int result(0);

    int lidx(p_localFromGlobalBusIndex(bus_idx));

    if (lidx >= 0) {
      result = p_network->getBus(lidx)->numStorage();
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
   * @param branch_idx @e global branch index
   * @return total number of lines in branch
   */
  int numLines(const int& branch_idx)
  {
    int result(0);

    int lidx(p_localFromGlobalBranchIndex(branch_idx));

    if (lidx >= 0) {
      result = p_network->getBranch(lidx)->numLines();
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
   * Collective
   *
   * @param bus_idx @e global bus index
   *
   * @result vector of @e global branch indexes
   **/
  std::vector<int> getConnectedBranches(const int& bus_idx) const
  {
    int root(0);
    std::vector<int> result;

    int lidx(p_localFromGlobalBusIndex(bus_idx));

    if (lidx >= 0) {
      typename NetworkType::BusPtr thebus(p_network->getBus(lidx));
      int oidx(thebus->getOriginalIndex());
      std::vector<int> lbranchi(p_network->getConnectedBranches(oidx));
      for (auto ibr : lbranchi) {
        result.push_back(p_network->getGlobalBranchIndex(ibr));
      }
      root = p_network->processor_rank();
    }

    p_network->communicator().sum(&root, 1);
    boost::mpi::broadcast(p_network->communicator(), result, root);

    return result;
  }

    
  /// Get the bus indexes connected to a branch
  /**
   * @e collective
   *
   * @param idx @e global branch index
   * @param fbus pointer to location to put from-bus (global) index
   * @param tbus pointer to location to put to-bus (global) index
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

  /**
   * Get a value from the bus' data collection.  To avoid coding this
   * a bunch of times, we'll use a broadcast, instead of
   * all_reduce, from the process where the active bus resides.
   */
  template <typename T>
  bool
  getBusInfo(const int& bus_idx, const std::string& field,
             T& value, const int& dev_idx = -1)
  {
    bool ok(false);
    int root(0);
    T result;

    int lidx(p_localFromGlobalBusIndex(bus_idx));

    if (lidx >= 0) {
      typename NetworkType::BusPtr thebus(p_network->getBus(lidx));
      if (dev_idx < 0) {
        ok = thebus->getData()->getValue(field.c_str(), &result);
      } else {
        ok = thebus->getData()->getValue(field.c_str(), &result, dev_idx);
      }
      if (ok) {
        root = p_network->processor_rank();
      }
    }

    ok = p_network->communicator().any(ok);
    if (ok) {
      p_network->communicator().sum(&root, 1);
      boost::mpi::broadcast(p_network->communicator(), result, root);
      value = result;
    }
    return ok;
  }

  /**
   * Get a value from the branch's data collection.  To avoid coding
   * this a bunch of times, we'll use a broadcast, instead of
   * all_reduce, from the process where the active branch resides.
   */
  template <typename T>
  bool
  getBranchInfo(const int& branch_idx, const std::string& field,
             T& value, const int& dev_idx = -1)
  {
    bool ok(false);
    int root(0);
    T result;

    int lidx(p_localFromGlobalBranchIndex(branch_idx));

    if (lidx >= 0) {
      typename NetworkType::BranchPtr thebranch(p_network->getBranch(lidx));
      if (dev_idx < 0) {
        ok = thebranch->getData()->getValue(field.c_str(), &result);
      } else {
        ok = thebranch->getData()->getValue(field.c_str(), &result, dev_idx);
      }
      if (ok) {
        root = p_network->processor_rank();
      }
    }

    ok = p_network->communicator().any(ok);
    if (ok) {
      p_network->communicator().sum(&root, 1);
      boost::mpi::broadcast(p_network->communicator(), result, root);
      value = result;
    }
    return ok;
  }
  
protected:

  /** 
   * Not collective.  
   *
   *Given a @e global bus index, find the local bus
   * index if it is active on this process.  
   * @param bus_idx global bus index
   * @return local bus index if active, -1 otherwise
   */ 
  int
  p_localFromGlobalBusIndex(const int& bus_idx) const
  {
    int result(-1);
    const int nBus(p_network->numBuses());
    for (int i = 0; i < nBus; ++i) {
      if (p_network->getBus(i)->getGlobalIndex() == bus_idx) {
        if (p_network->getActiveBus(i)) {
          result = i;
        }
      }
    }
    return result;
  }
  
  /** 
   * Not collective.  
   *
   * Given a @e global branch index, find the local branch
   * index if it is active on this process.  
   * @param bus_idx global branch index
   * @return local branch index if active, -1 otherwise
   */ 
  int
  p_localFromGlobalBranchIndex(const int& branch_idx) const
  {
    int result(-1);
    const int nBranch(p_network->numBranches());
    for (int i = 0; i < nBranch; ++i) {
      if (p_network->getBranch(i)->getGlobalIndex() == branch_idx) {
        if (p_network->getActiveBranch(i)) {
          result = i;
        }
      }
    }
    return result;
  }
  private:

  NetworkPtr p_network;
};

} // namespace gridpack
} // namespace analysis

#endif
