// Emacs Mode Line: -*- Mode:c++;-*-
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   global_vector.hpp
 * @author Bruce Palmer
 * @date   2017-07-5 09:25:03 d3g293
 * 
 * @brief  
 * This is a utility that is designed to store a collection of data elements
 * and make them accessible from any processor. The elements are assumed to
 * form a linear array.
 */

// -------------------------------------------------------------

#ifndef _global_vector_hpp_
#define _global_vector_hpp_

#include <iostream>
#include <ga.h>
#include "gridpack/utilities/exception.hpp"
#include "gridpack/parallel/communicator.hpp"

namespace gridpack {
namespace parallel {

// -------------------------------------------------------------
//  class GlobalVector
// -------------------------------------------------------------
template <typename _data_type >
 class GlobalVector {
private:
  gridpack::parallel::Communicator p_comm;
public:

  /**
   * Default constructor
   * @param comm communicator over which GlobalVector object runs.
   *             Data is accessible from any process on the communicator
   */
  GlobalVector(const gridpack::parallel::Communicator &comm)
    : p_comm(comm)
  {
    p_datasize = sizeof(_data_type);
    p_me = comm.rank();
    p_nprocs = comm.size();
    p_uploaded = false;
  }

  /**
   * Default destructor
   */
  ~GlobalVector(void)
  {
    if (p_uploaded) {
      NGA_Deregister_type(p_GA_Type);
      GA_Destroy(p_GA);
    }
  }

  /**
   * Add elements to GlobalVector
   * @param vec standard vector containing data
   * @param idx vector of indices
   */
  void addElements(const std::vector<int> &idx,
      const std::vector<_data_type> &vec)
  {
    int i, size;
    if (idx.size() != vec.size()) {
      printf("addElements: vector of indices does not match vector of data\n");
      return;
    }
    size = vec.size();
    for (i=0; i<size; i++) {
      p_index.push_back(idx[i]);
      p_data.push_back(vec[i]);
    }
  }

  /**
   *  Reset values of existing distributed array
   * @param vec standard vector containing data
   * @param idx vector of indices
   */
  void resetElements(const std::vector<int> &idx,
      const std::vector<_data_type> &vec)
  {
    int i, size;
    if (!p_uploaded) {
      printf("resetElements: Vector has not been previously uploaded!\n");
      return;
    }
    if (idx.size() != vec.size()) {
      printf("resetElements: vector of indices does not match vector of data\n");
      return;
    }
    size = vec.size();
    p_index.clear();
    p_data.clear();
    for (i=0; i<size; i++) {
      if (idx[i] >= p_numElems) {
        printf("resetElements: illegal index %d is greater than vector length %d\n",
            idx[i],p_numElems);
        continue;
      }
      p_index.push_back(idx[i]);
      p_data.push_back(vec[i]);
    }
  }

  /**
   * Upload data that is held locally into distributed array, so that it is
   * available anywhere in the system
   */
  void upload()
  {
    // Find out maximum index and assume that this represents the total number
    // of vectors added to GlobalVector object. If it doesn't, then there will be
    // some indices with no data.
    int max_idx = -1;
    int i;
    int one = 1;
    for (i=0; i<p_index.size(); i++) {
      if (p_index[i] < 0) {
        char buf[256];
        sprintf(buf,"Index %d less than zero in GlobalVector::upload on process %d\n",
            p_index[i],p_me);
        printf("%s",buf);
        throw gridpack::Exception(buf);
      }
      if (p_index[i] > max_idx) max_idx = p_index[i];
    }
    p_comm.max(&max_idx,1);
    p_numElems = max_idx+1;
    // Assume that if vector has no data then this is an error
    if (p_numElems < 1) {
      char buf[256];
      sprintf(buf,"Vector has no data in GlobalVector::upload on process %d\n",
          p_me);
      printf("%s",buf);
      throw gridpack::Exception(buf);
    }

    // Check to see if any indices have been used more than once. Throw an
    // exception if they have
    int g_idx = NGA_Create_handle();
    GA_Set_data(g_idx, one, &p_numElems, C_INT);
    int GAgrp = p_comm.getGroup();
    GA_Set_pgroup(g_idx, GAgrp);
    GA_Allocate(g_idx);
    int *iarr = new int[p_index.size()];
    int **iptr = new int*[p_index.size()];
    int *ival = new int[p_index.size()];
    GA_Zero(g_idx);
    for (i=0; i<p_index.size(); i++) {
      iarr[i] = p_index[i];
      iptr[i] = &iarr[i];
      ival[i] = one;
    }
    NGA_Scatter_acc(g_idx,ival,iptr,p_index.size(),&one);
    delete [] ival;
    int ilo, ihi, ld;
    int *lptr;
    NGA_Distribution(g_idx,p_me,&ilo,&ihi);
    NGA_Access(g_idx,&ilo,&ihi,&lptr,&ld);
    ld = ihi-ilo+1; 
    for (i=0; i<ld; i++) {
      if (lptr[i] > 1) {
        char buf[256];
        sprintf(buf,"GlobalVector::upload Multiple elements for index %d on process %d\n",
            i+ilo,p_me);
        printf("%s",buf);
        throw gridpack::Exception(buf);
      }
    }
    NGA_Release(g_idx,&ilo,&ihi);
    GA_Destroy(g_idx);

    // Create a GA to hold data
    p_GA_Type = NGA_Register_type(p_datasize);
    p_GA = GA_Create_handle();
    GA_Set_data(p_GA, one, &p_numElems, p_GA_Type);
    GA_Set_pgroup(p_GA, GAgrp);
    GA_Allocate(p_GA);

    // Make sure that vector is zeroed out
    _data_type *dptr;
    NGA_Distribution(p_GA,p_me,&ilo,&ihi);
    NGA_Access(p_GA,&ilo,&ihi,&dptr,&ld);
    ld = ihi-ilo+1;
    memset(dptr,'\0',ld*p_datasize);
    NGA_Release(p_GA,&ilo,&ihi);

    // Copy data to GA
    _data_type *data = new _data_type[p_index.size()];
    for (i=0; i<p_index.size(); i++) {
      data[i] = p_data[i];
    }
    NGA_Scatter(p_GA,data,iptr,p_index.size());
    delete [] iarr;
    delete [] iptr;
    delete [] data;
    p_data.clear();
    p_index.clear();
    p_uploaded = true;
    GA_Pgroup_sync(GAgrp);
  }

  /**
   * Upload new data values to existing distributed array.
   */
  void reload()
  {
    if (!p_uploaded) {
      printf("reload: Vector has not been previously uploaded!\n");
      return;
    }
    int i;
    int **iptr = new int*[p_index.size()];
    int *ival = new int[p_index.size()];
    for (i=0; i<p_index.size(); i++) {
      ival[i] = p_index[i];
      iptr[i] = &ival[i];
    }
    // Copy data to GA
    _data_type *data = new _data_type[p_index.size()];
    for (i=0; i<p_index.size(); i++) {
      data[i] = p_data[i];
    }
    NGA_Scatter(p_GA,data,iptr,p_index.size());
    delete [] ival;
    delete [] iptr;
    delete [] data;
    p_data.clear();
    p_index.clear();
    GA_Pgroup_sync(p_comm.getGroup());
  }

  /**
   * Get a vector vector of data elements corresponding to indices in idx
   * from the GlobalVector
   * @param idx vector of indices for requested data
   * @param vec vector of returned values
   */
  void getData(const std::vector<int> &idx, std::vector<_data_type> &vec)
  {
    int i;
    for (i=0; i<idx.size(); i++) {
      if (idx[i] < 0 || idx[i] >= p_numElems) {
        char buf[256];
        sprintf(buf,"GlobalVector::getData Requested vector index %d out of range on process %d\n",
            idx[i],p_me);
        printf("%s",buf);
        throw gridpack::Exception(buf);
      }
    }
    int *iarr = new int[idx.size()];
    int **iptr = new int*[idx.size()];
    for (i=0; i<idx.size(); i++) {
      iarr[i] = idx[i];
      iptr[i] = &iarr[i];
    }
    _data_type *dbuf = new _data_type[idx.size()];
    NGA_Gather(p_GA,dbuf,iptr,idx.size());
    // Copy results to vector
    vec.clear();
    for (int i=0; i<idx.size(); i++) {
      vec.push_back(dbuf[i]);
    }
    delete [] iarr;
    delete [] iptr;
    delete [] dbuf;
  };

  /**
   * Get all data from global vector
   * @param vec vector of returned values
   */
  void getAllData(std::vector<_data_type> &vec)
  {
    vec.clear();
    vec.resize(p_numElems);
    int lo, hi;
    int one = 1;
    lo = 0;
    hi = p_numElems-1;
    NGA_Get(p_GA,&lo,&hi,&vec[0],&one);
  }

private:
  // locally stored data before sending it to distributed array
  std::vector<_data_type> p_data;
  std::vector<int> p_index;

  // global array handle for storing distributed data
  int p_GA;

  // size of data type in bytes, processor configuration information,
  // total number stored vectors
  int p_datasize;
  int p_me;
  int p_nprocs;
  int p_numElems;
  int p_GA_Type;

  // flag to track if data has been uploaded
  bool p_uploaded;
};
} // namespace parallel
} // namespace gridpack

#endif

