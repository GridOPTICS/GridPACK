// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hash_distr.hpp
 * @author Bruce Palmer
 * @date   2014-07-28 09:25:03 d3g293
 * 
 * @brief  
 * This is a utility that is designed to move a set of data for a collection of
 * buses and/or branches based on the process that owns the buses or branches.
 * Ownership is based on the original index of the bus or branch
 * 
 */

// -------------------------------------------------------------

#ifndef _hash_distr_hpp_
#define _hash_distr_hpp_

#include "gridpack/parallel/index_hash.hpp"

namespace gridpack {
namespace hash_distr {

// -------------------------------------------------------------
//  class HashDistribution
// -------------------------------------------------------------
template <typename _network,
          typename _bus_data_type,
          typename _branch_data_type>
 class HashDistribution {
private:

  typedef struct{bool flag;
           _bus_data_type data;
  } bus_data_pair;

  typedef struct{bool flag;
           _branch_data_type data;
  } branch_data_pair;

public:
  typedef _network NetworkType;
  typedef boost::shared_ptr<NetworkType> NetworkPtr;

  // Default constructor
  // @param network network over which hash distribution function extends
  HashDistribution(const boost::shared_ptr<_network> network)
    : p_network(network)
  {
    p_hashMap.reset(new
        gridpack::hash_map::GlobalIndexHashMap(p_network->communicator()));
    p_GAgrp = p_network->communicator().getGroup();

    // Initialize hash map using original bus indices and global indices
    int i,ikey,ival,idx1,idx2;
    std::vector<std::pair<int,int> > busPairs;
    int nbus = p_network->numBuses();
    for (i=0; i<nbus; i++) {
      ikey = p_network->getOriginalBusIndex(i);
      ival = p_network->getGlobalBusIndex(i);
      busPairs.push_back(std::pair<int,int>(ikey,ival));
    }
    p_hashMap->addPairs(busPairs);
    busPairs.clear();
    std::pair<int,int> key;
    std::vector<std::pair<std::pair<int,int>, int> >
      branchPairs;
    int nbranch = p_network->numBranches();
    for (i=0; i<nbranch; i++) {
      p_network->getOriginalBranchEndpoints(i,&idx1,&idx2);
      key = std::pair<int,int>(idx1,idx2);
      ival =
        p_network->getGlobalBranchIndex(i);
      branchPairs.push_back(std::pair<std::pair<int,int>,int>(key,ival));
    }
    p_hashMap->addPairs(branchPairs);
    branchPairs.clear();

    //store size of data items
    p_size_bus_data = sizeof(bus_data_pair);
    p_size_branch_data = sizeof(branch_data_pair);
  }

  // Default destructor
  ~HashDistribution(void)
  {
  }

  // Send values corresponding to keys to the processors that own them. On
  // completion, keys contain the local indices of the buses recieving data and
  // values contains the data
  // @param keys on input, a list of integer keys corresponding to the original
  // bus indices of the buses that receive the data, on output, a list of local
  // bus indices
  // @param values on input, a list of values corresponding to the list of keys,
  // on output, the values of the received data
  void distributeBusValues(std::vector<int> &keys, std::vector<_bus_data_type> &values)
  {
    // Get global indices corresponding to keys
    std::vector<int> g_idx;
    p_hashMap->getValues(keys,g_idx);

    // Create global array to distribute data
    int g_bus = GA_Create_handle();
    int g_type = NGA_Register_type(p_size_bus_data);
    int one = 1;
    int totalBuses = p_network->totalBuses();
    GA_Set_data(g_bus,one,&totalBuses,g_type);
    GA_Set_pgroup(g_bus,p_GAgrp);
    GA_Allocate(g_bus);
    NGA_Deregister_type(g_type);

    // initialize all data pairs in global array to false
    int i, nsize, lo, hi, ld;
    int me = GA_Nodeid();
    bus_data_pair *list;
    NGA_Distribution(g_bus,me,&lo,&hi);
    nsize = hi - lo + 1;
    NGA_Access(g_bus,&lo,&hi,&list,&ld);
    for (i=0; i<nsize; i++) {
      list[i].flag = false;
    }
    NGA_Release(g_bus,&lo,&hi);

    // Copy indices and values to local arrays
    nsize = values.size();
    int *index = new int[nsize];
    int **idx = new int*[nsize];
    bus_data_pair *values_buf = new bus_data_pair[nsize];
    for (i=0; i<nsize; i++) {
      index[i] = g_idx[i];
      idx[i] = &index[i];
      values_buf[i].flag = true;
      values_buf[i].data = values[i];
    }
    // Scatter data to global array
    NGA_Scatter(g_bus,values_buf,idx,nsize);
    GA_Pgroup_sync(p_GAgrp);
    delete [] index;
    delete [] idx;
    delete [] values_buf;
    keys.clear();
    values.clear();

    // Set up local arrays to receive data
    nsize = p_network->numBuses();
    index = new int[nsize];
    idx = new int*[nsize];
    values_buf = new bus_data_pair[nsize];
    int icnt = 0;
    for (i=0; i<nsize; i++) {
      if (p_network->getActiveBus(i)) {
        index[icnt] = p_network->getGlobalBusIndex(i);
        idx[icnt] = &index[icnt];
        icnt++;
      }
    }

    // Gather data to local buffers
    NGA_Gather(g_bus,values_buf,idx,icnt);

    // Copy data back to vectors
    icnt = 0;
    for (i=0; i<nsize; i++) {
      if (p_network->getActiveBus(i)) {
        if (values_buf[icnt].flag) {
          keys.push_back(i);
          values.push_back(values_buf[icnt].data);
        }
        icnt++;
      }
    }
    delete [] index;
    delete [] idx;
    delete [] values_buf;
    GA_Destroy(g_bus);
  }

  // Send values corresponding to keys to the processors that own them. On
  // completion, keys contain the local indices of the branches recieving data
  // and
  // values contains the data
  // @param keys on input, a list of integer pair keys corresponding to the
  // branches that receive the data, on output, a list of local branch indices
  // @param values on input, a list of values corresponding to the list of keys,
  // on output, the values of the received data
  void distributeBranchValues(std::vector<std::pair<int,int> > &keys,
      std::vector<int> &branch_ids,
      std::vector<_branch_data_type> &values)
  {
    // Get global indices corresponding to keys
    std::vector<int> g_idx;
    p_hashMap->getValues(keys,g_idx);
    // Create global array to distribute data
    int g_branch = GA_Create_handle();
    int g_type = NGA_Register_type(p_size_branch_data);
    int one = 1;
    int totalBranches = p_network->totalBranches();
    GA_Set_data(g_branch,one,&totalBranches,g_type);
    GA_Set_pgroup(g_branch,p_GAgrp);
    GA_Allocate(g_branch);
    NGA_Deregister_type(g_type);

    // initialize all data pairs in global array to false
    int i, nsize, lo, hi, ld;
    int me = GA_Nodeid();
    branch_data_pair *list;
    NGA_Distribution(g_branch,me,&lo,&hi);
    NGA_Access(g_branch,&lo,&hi,&list,&ld);
    nsize = hi - lo + 1;
    for (i=0; i<nsize; i++) {
      list[i].flag = false;
    }
    NGA_Release(g_branch,&lo,&hi);

    // Copy keys and values to local arrays
    nsize = values.size();
    int *index = new int[nsize];
    int **idx = new int*[nsize];
    branch_data_pair *values_buf = new branch_data_pair[nsize];
    for (i=0; i<nsize; i++) {
      index[i] = g_idx[i];
      idx[i] = &index[i];
      values_buf[i].flag = true;
      values_buf[i].data = values[i];
    }
    // Scatter data to global array
    NGA_Scatter(g_branch,values_buf,idx,nsize);
    GA_Pgroup_sync(p_GAgrp);
    delete [] index;
    delete [] idx;
    delete [] values_buf;
    g_idx.clear();
    values.clear();

    // Set up local arrays to receive data
    nsize = p_network->numBranches();
    index = new int[nsize];
    idx = new int*[nsize];
    values_buf = new branch_data_pair[nsize];
    int icnt = 0;
    for (i=0; i<nsize; i++) {
      if (p_network->getActiveBranch(i)) {
        index[icnt] = p_network->getGlobalBranchIndex(i);
        idx[icnt] = &index[icnt];
        icnt++;
      }
    }

    // Gather data to local buffers
    NGA_Gather(g_branch,values_buf,idx,icnt);

    // Copy data back to vectors
    icnt = 0;
    for (i=0; i<nsize; i++) {
      if (p_network->getActiveBranch(i)) {
        if (values_buf[icnt].flag) {
          branch_ids.push_back(i);
          values.push_back(values_buf[icnt].data);
        }
        icnt++;
      }
    }
    delete [] index;
    delete [] idx;
    delete [] values_buf;
    GA_Destroy(g_branch);
  }

private:

  boost::shared_ptr<gridpack::hash_map::GlobalIndexHashMap> p_hashMap;

  NetworkPtr p_network;

  int p_size_bus_data;
  int p_size_branch_data;

  int p_GAgrp;

};


} // namespace gridpack
} // namespace hash_distr

#endif

