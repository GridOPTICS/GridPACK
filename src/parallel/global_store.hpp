// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   global_store.hpp
 * @author Bruce Palmer
 * @date   2016-07-8 09:25:03 d3g293
 * 
 * @brief  
 * This is a utility that is designed to store a collection of vectors and make
 * them accessible from any processor
 */

// -------------------------------------------------------------

#ifndef _global_store_hpp_
#define _global_store_hpp_

#include <iostream>
#include <ga.h>
#include "gridpack/utilities/exception.hpp"
#include "gridpack/parallel/communicator.hpp"

namespace gridpack {
namespace parallel {

// -------------------------------------------------------------
//  class GlobalStore
// -------------------------------------------------------------
template <typename _data_type >
 class GlobalStore {
private:
  gridpack::parallel::Communicator p_comm;
public:

  /**
   * Default constructor
   * @param comm communicator over which GlobalStore object runs.
   *             Data is accessible from any process on the communicator
   */
  GlobalStore(const gridpack::parallel::Communicator &comm)
    : p_comm(comm)
  {
    p_begin = NULL;
    p_end = NULL;
    p_datasize = sizeof(_data_type);
    p_me = comm.rank();
    p_nprocs = comm.size();
    p_uploaded = false;
  }

  /**
   * Default destructor
   */
  ~GlobalStore(void)
  {
    // Assume that if p_begin is allocated, other data objects are allocated
    if (p_uploaded) {
      GA_Destroy(p_GA);
    }
    if (p_begin != NULL) delete [] p_begin;
    if (p_end != NULL) delete [] p_end;
  }

  /**
   * Add vector to GlobalStore
   * @param vec standard vector containing data
   * @param idx index of data in GlobalStore
   */
  void addVector(const int idx, const std::vector<_data_type> &vec)
  {
    p_index.push_back(idx);
    p_data.push_back(vec);
  }

  /**
   * Upload data that is held locally into distributed array, so that it is
   * available anywhere in the system
   */
  void upload()
  {
    // Find out maximum index and assume that this represents the total number
    // of vectors added to GlobalStore object. If it doesn't, then there will be
    // some indices with no data.
    int max_idx = 0;
    int i;
    for (i=0; i<p_index.size(); i++) {
      if (p_index[i] < 0) {
        char buf[256];
        sprintf(buf,"Index %d less than zero in GlobalStore::upload on process %d\n",
            p_index[i],p_me);
        printf("%s",buf);
        throw gridpack::Exception(buf);
      }
      if (p_index[i] > max_idx) max_idx = p_index[i];
    }
    p_comm.max(&max_idx,1);
    p_numVecs = max_idx+1;

    // Allocate memory for index arrays marking start and end of data
    int nsize = max_idx+1;
    p_begin = new int[nsize];
    p_end = new int[nsize];
    for (i=0; i<nsize; i++) {
      p_begin[i] = 0;
      p_end[i] = 0;
    }

    // Check to see if any indices have been used more than once. Throw an
    // exception if they have
    for (i=0; i<p_index.size(); i++) {
      p_begin[p_index[i]] += 1;
    }
    p_comm.sum(p_begin,nsize);
    for (i=0; i<nsize; i++) {
      if (p_begin[i] > 1) {
        char buf[256];
        sprintf(buf,"Multiple vectors for index %d GlobalStore::upload on process %d\n",
            i,p_me);
        printf("%s",buf);
        throw gridpack::Exception(buf);
      }
      p_begin[i] = 0;
    }
    //  Broadcast size of individual vectors to all processors. These can be
    //  used to construct the begin and end arrays as well determine the size of
    //  the global array used for distributed storage
    for (i=0; i<p_index.size(); i++) {
      p_end[p_index[i]] = p_data[i].size();
    }
    p_comm.sum(p_end,nsize);
    p_begin[0] = 0;
    for (i=1; i<nsize; i++) {
      p_begin[i] = p_begin[i-1]+p_end[i-1];
    }
    int ndata = 0;
    for (i=0; i<nsize; i++) {
      ndata += p_end[i];
      p_end[i] = p_begin[i]+p_end[i];
    }

    if (ndata < 1) {
      std::cout << "No data found for global store" << std::endl;
      return;
    }

    // Create a GA to hold data
    int one = 1;
    int GA_type = NGA_Register_type(p_datasize);
    p_GA = GA_Create_handle();
    GA_Set_data(p_GA, one, &ndata, GA_type);
    int GAgrp = p_comm.getGroup();
    GA_Set_pgroup(p_GA, GAgrp);
    GA_Allocate(p_GA);

    // Copy data to GA
    for (i=0; i<p_index.size(); i++) {
      int idx = p_index[i];
      int lo = p_begin[idx];
      int hi = p_end[idx]-1;
      NGA_Put(p_GA,&lo,&hi,&((p_data[i])[0]),&one);
      p_data[i].clear();
    }
    p_data.clear();
    p_index.clear();
    p_uploaded = true;
    GA_Pgroup_sync(GAgrp);
  }

  /**
   * Get vector corresponding to index idx from GlobalStore
   * @param idx index of stored vector
   * @param vec vector of returned values
   */
  void getVector(const int idx, std::vector<_data_type> &vec)
  {
    if (idx < 0 || idx >= p_numVecs) {
      char buf[256];
      sprintf(buf,"Requested vector index %d in GlobalStore::getVector out of range on process %d\n",
          idx,p_me);
      printf("%s",buf);
      throw gridpack::Exception(buf);
    }
    int ld = 1;
    int lo = p_begin[idx];
    int hi = p_end[idx]-1;
    int size = p_end[idx]-p_begin[idx];
    _data_type *dbuf = new _data_type[size];
    NGA_Get(p_GA,&lo,&hi,dbuf,&ld);
    // Copy results to vector
    vec.clear();
    for (int i=0; i<size; i++) {
      vec.push_back(dbuf[i]);
    }
    delete [] dbuf;
  };

private:
  // beginning and end indices of stored data
  int *p_begin;
  int *p_end;

  // locally stored data before sending it to distributed array
  std::vector<std::vector<_data_type> > p_data;
  std::vector<int> p_index;

  // global array handle for storing distributed data
  int p_GA;

  // size of data type in bytes, processor configuration information,
  // total number stored vectors
  int p_datasize;
  int p_me;
  int p_nprocs;
  int p_numVecs;

  // flag to track if data has been uploaded
  bool p_uploaded;
};
} // namespace gridpack
} // namespace utility

#endif

