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

#define SYSTOLIC

//#define HASH_WITH_MPI

#include <ga.h>
#include <boost/unordered_map.hpp>
#include "gridpack/parallel/index_hash.hpp"
#include "gridpack/utilities/exception.hpp"

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

  typedef struct{int idx;
           _bus_data_type data;
  } bus_data_pair;

  typedef struct{int idx1;
            int idx2;
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
#ifndef SYSTOLIC
    p_indexHashMap.reset(new
        gridpack::hash_map::GlobalIndexHashMap(p_network->communicator()));
#endif
    p_GAgrp = p_network->communicator().getGroup();
    int me = p_network->communicator().rank();

#ifndef SYSTOLIC
    // Initialize hash map using original bus indices and global indices
    int i,ikey,ival,idx1,idx2;
    std::vector<std::pair<int,int> > busPairs;
    int nbus = p_network->numBuses();
    for (i=0; i<nbus; i++) {
      ikey = p_network->getOriginalBusIndex(i);
      ival = me;
      busPairs.push_back(std::pair<int,int>(ikey,ival));
    }
    p_indexHashMap->addPairs(busPairs);
    busPairs.clear();
    std::pair<int,int> key;
    std::vector<std::pair<std::pair<int,int>, int> >
      branchPairs;
    int nbranch = p_network->numBranches();
    for (i=0; i<nbranch; i++) {
      p_network->getOriginalBranchEndpoints(i,&idx1,&idx2);
      key = std::pair<int,int>(idx1,idx2);
      ival = me;
      branchPairs.push_back(std::pair<std::pair<int,int>,int>(key,ival));
    }
    p_indexHashMap->addPairs(branchPairs);
    branchPairs.clear();
#endif

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
#ifdef SYSTOLIC
    int ksize = keys.size();
    int vsize = values.size();
    int me = GA_Pgroup_nodeid(p_GAgrp);
    int nprocs = GA_Pgroup_nnodes(p_GAgrp);
    if (vsize != ksize) {
      char buf[256];
      sprintf(buf,"p[%d] HashDistribution::distributeBusValues ERROR: length of"
          " keys and values arrays don't match ksize: %d vsize: %d\n",
          me,ksize,vsize);
      printf("%s",buf);
      throw gridpack::Exception(buf);
    }
    
    // Construct mapc array for irregular distribution
    int i;
    int *sizes = new int[nprocs];
    for (i=0; i<nprocs; i++) {
      sizes[i] = 0;
    }
    sizes[me] = ksize;
    char plus[2];
    strcpy(plus,"+");
    GA_Pgroup_igop(p_GAgrp,sizes,nprocs,plus);
    int *mapc = new int[nprocs];
    mapc[0] = 0;
    int total_values = sizes[0];
    for (i=1; i<nprocs; i++) {
      mapc[i] = mapc[i-1]+sizes[i-1];
      total_values += sizes[i];
    }
    if (total_values == 0) {
      delete [] sizes;
      delete [] mapc;
      return;
    }

    // Copy values to a straight array
    bus_data_pair *list;
    if (ksize > 0) {
      list = new bus_data_pair[ksize];
      for (i=0; i<ksize; i++) {
        list[i].idx = keys[i];
        list[i].data = values[i];
      }
    }
    int g_type = NGA_Register_type(p_size_bus_data);

    // Create global array and store all values in it
    int lo, hi;
    lo = 0;
    for (i=0; i<me; i++) {
      lo += sizes[i];
    }
    hi = lo + sizes[me] - 1;
    int one = 1;
    int g_vals = GA_Create_handle();
    GA_Set_data(g_vals,one,&total_values,g_type);
//    GA_Set_irreg_distr(g_vals,mapc,&nprocs);
    GA_Set_pgroup(g_vals,p_GAgrp);
    if (!GA_Allocate(g_vals)) {
      char buf[256];
      sprintf(buf,"HashDistribution::distributeBusValues: Unable to allocate"
          " distributed array for storing values");
      printf("%s",buf);
      throw gridpack::Exception(buf);
    }
    if (lo <= hi) NGA_Put(g_vals, &lo, &hi, list, &one);
    GA_Pgroup_sync(p_GAgrp);
    NGA_Deregister_type(g_type);
    if (ksize > 0) delete [] list;
    delete [] mapc;
    delete [] sizes;

    // Create a map that maps original bus indices to local indices
    int idx;
    int nbus = p_network->numBuses();
    std::multimap<int,int> hmap;
    for (i=0; i<nbus; i++) {
      idx = p_network->getOriginalBusIndex(i);
      hmap.insert(std::pair<int,int>(idx,i));
    }
    std::multimap<int,int>::iterator it;
    
    // Get values from global array and copy them into output arrays
    keys.clear();
    values.clear();
    double delta = static_cast<double>(total_values)/static_cast<double>(nprocs);
    for (i = 0; i<nprocs; i++) {
      idx = (i+me)%nprocs;
      lo = static_cast<int>(delta*static_cast<double>(idx));
      if (idx<nprocs-1) {
        hi = static_cast<int>(delta*static_cast<double>(idx+1))-1;
      } else {
        hi = total_values-1;
      }
      int nsize = hi - lo + 1;
      if (lo <= hi) {
        list = new bus_data_pair[nsize];
        if (lo<=hi) NGA_Get(g_vals, &lo, &hi, list, &one);
        int j;
        for (j=0; j<nsize; j++) {
          it = hmap.find(list[j].idx);
          if (it != hmap.end()) {
            while (it != hmap.upper_bound(list[j].idx)) {
              keys.push_back(it->second);
              values.push_back(list[j].data);
              it++;
            }
          }
        }
        delete [] list;
      }
    }
    GA_Destroy(g_vals);
#else
    int nprocs = p_network->communicator().size();
    int me = p_network->communicator().rank();
    // Get global indices corresponding to keys
    int i,j;
    // Create a list of unique keys
    std::vector<int> base_keys;
    std::set<int> key_check;
    std::set<int>::iterator itc;
    for (i=0; i<keys.size(); i++) {
      itc = key_check.find(keys[i]);
      if (itc == key_check.end()) {
        key_check.insert(keys[i]);
        base_keys.push_back(keys[i]);
      }
    }
    std::vector<int> procLoc;
    // Get processor that owns global index
    p_indexHashMap->getValues(base_keys,procLoc);
    // We now know how many processors each unique key maps to. Use this
    // information to create a complete a  new values list, a new keys list, and
    // a processor location list
    std::multimap<int,int> keyMap;
    for (i=0; i<base_keys.size(); i++) {
      keyMap.insert(std::pair<int,int>(base_keys[i],procLoc[i]));
    }
    std::vector<_bus_data_type> newValues;
    std::vector<int> newKeys;
    std::vector<int> destProcs;
    j = 0;
    std::multimap<int,int>::iterator itk;
    for (i=0; i<values.size(); i++) {
      itk = keyMap.find(keys[i]);
      if (itk != keyMap.end()) {
        while (itk != keyMap.upper_bound(keys[i])) {
          newValues.push_back(values[i]);
          newKeys.push_back(keys[i]);
          destProcs.push_back(itk->second);
          itk++;
        }
      }
    }

    // sort keys into linked list based on which processors they are going to
    int nsize = newKeys.size();
    int ltop[nprocs];
    int ldest[nsize];
    int destNum[nprocs];
    for (i=0; i<nprocs; i++) {
      ltop[i] = -1;
      destNum[i] = 0;
    }
    for (i=0; i<nsize; i++) {
      ldest[i] = -1;
    }
    int iproc;
    for (i=0; i<nsize; i++) {
      iproc = destProcs[i];
      destNum[iproc]++;
      ldest[i] = ltop[iproc]; 
      ltop[iproc] = i;
    }
#ifdef HASH_WITH_MPI
    // Use all-to-all call to determine number of unique keys coming from each
    // processor
    int srcNum[nprocs];
    int ierr;
    int one = 1;
    MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
    ierr = MPI_Alltoall(destNum,one,MPI_INT,srcNum,one,MPI_INT,comm);

    // Each process now knows how much data it will receive. Pack data data into
    // an appropriate sized buffer and send it to processors using a all-to-all
    // vector call
    bus_data_pair *sendBuf;
    sendBuf = new bus_data_pair[newValues.size()];
    // Fill sendBuf with data from newValues and newKeys
    int icnt=0;
    for (i=0; i<nprocs; i++) {
      j = ltop[i];
      if (j>=0) {
        while(j >= 0) {
          sendBuf[icnt].idx = newKeys[j];
          sendBuf[icnt].data = newValues[j];
          j = ldest[j];
          icnt++;
        }
      }
    }

    // Evaluate offsets for data being sent and received and then correct all
    // sizes by the number of bytes in bus_data_pair structure
    int elemsize =  sizeof(bus_data_pair);
    int destOffset[nprocs];
    int srcOffset[nprocs];
    destOffset[0] = 0;
    srcOffset[0] = 0;
    for (i=1; i<nprocs; i++) {
      destOffset[i] = destOffset[i-1] + destNum[i-1];
      srcOffset[i] = srcOffset[i-1] + srcNum[i-1];
    }
    int nvalues = 0;
    for (i=0; i<nprocs; i++) {
      nvalues += srcNum[i];
      destOffset[i] = elemsize*destOffset[i];
      destNum[i] = elemsize*destNum[i];
      srcOffset[i] = elemsize*srcOffset[i];
      srcNum[i] = elemsize*srcNum[i];
    }

    // Allocate buffer to receive data from other processors
    bus_data_pair *recvBuf;
    recvBuf = new bus_data_pair[nvalues];
    
    // Transmit data and clean up buffers that are no longer needed
    ierr = MPI_Alltoallv(sendBuf, destNum, destOffset, MPI_BYTE, recvBuf,
        srcNum, srcOffset, MPI_BYTE, comm);
    delete [] sendBuf;

    // Data is now available on processor that can use it. Pack it into final
    // data structures. Start by creating a map between global and local bus
    // indices on this processors
    keys.clear();
    values.clear();
    int nbus = p_network->numBuses();
    std::multimap<int,int> idxMap;
    std::pair<int,int> idxPair;
    // Create map between global index and local index for locally held buses
    for (i=0; i<nbus; i++) {
      idxPair = std::pair<int,int>(p_network->getOriginalBusIndex(i),i);
      idxMap.insert(idxPair);
    }
    // Get data from recvBuf and pack keys and values arrays
    // with local indices and data
    std::multimap<int,int>::iterator it;
    for (i=0; i<nvalues; i++) {
      it = idxMap.find(recvBuf[i].idx);
      if (it != idxMap.end()) {
        while (it != idxMap.upper_bound(recvBuf[i].idx)) {
          keys.push_back(it->second);
          values.push_back(recvBuf[i].data);
          it++;
        }
      } else {
        printf("p[%d] Unresolved original bus index: %d\n",me,
            recvBuf[i].idx);
      }
    }

    delete [] recvBuf;
#else
    // Find total number of values going to each processor
    // and get offsets on each processor for each set  of values
    int g_numValues = GA_Create_handle();
    int dims;
    dims = nprocs;
    int one = 1;
    int blocks;
    blocks = 1;
    GA_Set_data(g_numValues,one,&dims,C_INT);
    GA_Set_chunk(g_numValues, &blocks);
    GA_Set_pgroup(g_numValues, p_GAgrp);
    GA_Allocate(g_numValues);
    GA_Zero(g_numValues);
    int r_offset[nprocs];
    for (j=0; j<nprocs; j++) {
      i = (j+me)%nprocs;
      if (destNum[i] > 0) {
        r_offset[i] = NGA_Read_inc(g_numValues,&i,destNum[i]);
      } else {
        r_offset[i] = 0;
      }
    }
    GA_Pgroup_sync(p_GAgrp);
    // Copy values from global array to local array
    int numValues[nprocs];
    int lo, hi;
    lo = 0;
    hi = nprocs-1;
    if (me == 0) {
      if (lo<=hi) NGA_Get(g_numValues,&lo,&hi,numValues,&one);
    } else {
      for (i=0; i<nprocs; i++) {
        numValues[i] = 0;
      }
    }
    char plus[2];
    strcpy(plus,"+");
    GA_Pgroup_igop(p_GAgrp,numValues,nprocs,plus);
    GA_Destroy(g_numValues);

    // Create a global array that can hold all values. Partition the array so
    // that each processor has enough local memory to hold all values coming to
    // it
    int dtype = NGA_Register_type(sizeof(bus_data_pair));
    // total number of values on all processors
    int totalVals = 0;
    for (i=0; i<nprocs; i++) {
      r_offset[i] += totalVals;
      totalVals += numValues[i];
    }
    if (totalVals == 0) {
      NGA_Deregister_type(dtype);
      return;
    }
    int g_data = GA_Create_handle();
    dims = totalVals;
    GA_Set_data(g_data, one, &dims, dtype);
    GA_Set_pgroup(g_data,p_GAgrp);
    int mapc[nprocs];
    mapc[0] = 0;
    for (i=1; i<nprocs; i++) {
      mapc[i] = mapc[i-1] + numValues[i-1];
    }
    blocks = nprocs;
    GA_Set_irreg_distr(g_data,mapc,&blocks);
    GA_Allocate(g_data);

    // Repack values and send them to the processor with the corresponding buses
    bus_data_pair *bus_data;
    int ncnt;
    for (i=0; i<nprocs; i++) {
      j = ltop[i];
      ncnt = 0;
      if (j >= 0) {
        bus_data = new bus_data_pair[destNum[i]];
        while (j >= 0) {
          bus_data[ncnt].idx = newKeys[j];
          bus_data[ncnt].data = newValues[j];
          j = ldest[j];
          ncnt++;
        }
        lo = r_offset[i];
        hi = lo + destNum[i] - 1;
        if (lo<=hi) NGA_Put(g_data,&lo,&hi,bus_data,&one);
        delete [] bus_data;
      }
    }

    GA_Pgroup_sync(p_GAgrp);
    // Data is now on processor. Unpack it and put it in arrays for export
    keys.clear();
    values.clear();
    int nbus = p_network->numBuses();
    std::multimap<int,int> idxMap;
    std::pair<int,int> idxPair;
    // Create map between global index and local index for locally held buses
    for (i=0; i<nbus; i++) {
      idxPair = std::pair<int,int>(p_network->getOriginalBusIndex(i),i);
      idxMap.insert(idxPair);
    }
    // Get data from global array and pack keys and values arrays
    // with local indices and data
    int ndata = numValues[me];
    lo = mapc[me];
    hi = lo + ndata - 1;
    if (lo<=hi) NGA_Access(g_data,&lo,&hi,&bus_data,&one);
    std::multimap<int,int>::iterator it;
    for (i=0; i<ndata; i++) {
      it = idxMap.find(bus_data[i].idx);
      if (it != idxMap.end()) {
        while (it != idxMap.upper_bound(bus_data[i].idx)) {
          keys.push_back(it->second);
          values.push_back(bus_data[i].data);
          it++;
        }
      } else {
        printf("p[%d] Unresolved original bus index: %d\n",me,
            bus_data[i].idx);
      }
    }
    if (lo<=hi) NGA_Release(g_data,&lo,&hi);
    GA_Destroy(g_data);
#endif
#endif
  }

  // Send values corresponding to keys to the processors that own them. On
  // completion, keys contain the local indices of the buses recieving data and
  // values contains the data
  // @param keys on input, a list of integer keys corresponding to the original
  // bus indices of the buses that receive the data, on output, a list of local
  // bus indices
  // @param values on input, a vector of pointers to arrays of values. Each
  // array is associeted with a bus. On output, the arrays ov values of the received data
  // @param nvals the number of values in each array
  void distributeBusValues(std::vector<int> &keys, std::vector<_bus_data_type*>
      &values, int nvals)
  {
#ifdef SYSTOLIC
    int ksize = keys.size();
    int vsize = values.size();
    int me = GA_Pgroup_nodeid(p_GAgrp);
    int nprocs = GA_Pgroup_nnodes(p_GAgrp);
    if (vsize != ksize) {
      char buf[256];
      sprintf(buf,"p[%d] HashDistribution::distributeBusValues ERROR: length"
          " of keys and values arrays don't match ksize: %d vsize: %d\n",
          me,ksize,vsize);
      printf("%s",buf);
      throw gridpack::Exception(buf);
    }
    
    // Construct mapc array for irregular distribution
    int i, j;
    int *sizes = new int[nprocs];
    for (i=0; i<nprocs; i++) {
      sizes[i] = 0;
    }
    sizes[me] = ksize;
    char plus[2];
    strcpy(plus,"+");
    GA_Pgroup_igop(p_GAgrp,sizes,nprocs,plus);
    int *mapc = new int[nprocs];
    mapc[0] = 0;
    int total_values = sizes[0];
    for (i=1; i<nprocs; i++) {
      mapc[i] = mapc[i-1]+sizes[i-1];
      total_values += sizes[i];
    }
    if (total_values == 0) {
      delete [] sizes;
      delete [] mapc;
      return;
    }

    // Copy values to a straight array
    char *list, *ptr;
    if (ksize > 0) {
      list = new char[ksize*(nvals*sizeof(_bus_data_type)+sizeof(int))];
      ptr = list;
      for (i=0; i<ksize; i++) {
        ((int*)ptr)[0] = keys[i];
        ptr += sizeof(int);
        for (j=0; j<nvals; j++) {
          ((_bus_data_type*)ptr)[j] = (values[i])[j];
        }
        ptr += nvals*sizeof(_bus_data_type);
      }
    }
    int g_type = NGA_Register_type(nvals*sizeof(_bus_data_type)+sizeof(int));

    // Create global array and store all values in it
    int lo, hi;
    lo = 0;
    for (i=0; i<me; i++) {
      lo += sizes[i];
    }
    hi = lo + sizes[me] - 1;
    int one = 1;
    int g_vals = GA_Create_handle();
    GA_Set_data(g_vals,one,&total_values,g_type);
//    GA_Set_irreg_distr(g_vals,mapc,&nprocs);
    GA_Set_pgroup(g_vals,p_GAgrp);
    if (!GA_Allocate(g_vals)) {
      char buf[256];
      sprintf(buf,"HashDistribution::distributeBusValues: Unable to allocate"
          " distributed array for storing values");
      printf("%s",buf);
      throw gridpack::Exception(buf);
    }
    if (lo <= hi) NGA_Put(g_vals, &lo, &hi, list, &one);
    GA_Pgroup_sync(p_GAgrp);
    NGA_Deregister_type(g_type);
    if (ksize > 0) delete [] list;
    delete [] mapc;
    delete [] sizes;

    // Create a map that maps original bus indices to local indices
    int idx;
    int nbus = p_network->numBuses();
    std::multimap<int,int> hmap;
    for (i=0; i<nbus; i++) {
      idx = p_network->getOriginalBusIndex(i);
      hmap.insert(std::pair<int,int>(idx,i));
    }
    std::multimap<int,int>::iterator it;
    
    // Get values from global array and copy them into output arrays
    keys.clear();
    for (i=0; i<ksize; i++) {
      delete [] values[i];
    }
    values.clear();
    double delta = static_cast<double>(total_values)/static_cast<double>(nprocs);
    for (i = 0; i<nprocs; i++) {
      idx = (i+me)%nprocs;
      lo = static_cast<int>(delta*static_cast<double>(idx));
      if (idx<nprocs-1) {
        hi = static_cast<int>(delta*static_cast<double>(idx+1))-1;
      } else {
        hi = total_values-1;
      }
      int nsize = hi - lo + 1;
      if (lo <= hi) {
        list = new char[nsize*(nvals*sizeof(_bus_data_type)+sizeof(int))];
        if (lo<=hi) NGA_Get(g_vals, &lo, &hi, list, &one);
        int j, k;
        ptr = list;
        for (j=0; j<nsize; j++) {
          idx = ((int*)ptr)[0];
          ptr += sizeof(int);
          it = hmap.find(idx);
          if (it != hmap.end()) {
            while (it != hmap.upper_bound(idx)) {
              keys.push_back(it->second);
              _bus_data_type *data = new _bus_data_type[nvals];
              for (k=0; k<nvals; k++) {
                data[k] = ((_bus_data_type*)ptr)[k];
              }
              values.push_back(data);
              it++;
            }
          }
          ptr += nvals*sizeof(_bus_data_type);
        }
        delete [] list;
      }
    }
    GA_Destroy(g_vals);
#else
    int nprocs = p_network->communicator().size();
    int me = p_network->communicator().rank();
    // Get global indices corresponding to keys
    int i,j;
    // Create a list of unique keys
    std::vector<int> base_keys;
    std::set<int> key_check;
    std::set<int>::iterator itc;
    for (i=0; i<keys.size(); i++) {
      itc = key_check.find(keys[i]);
      if (itc == key_check.end()) {
        key_check.insert(keys[i]);
        base_keys.push_back(keys[i]);
      }
    }
    std::vector<int> procLoc;
    // Get processor that owns global index
    p_indexHashMap->getValues(base_keys,procLoc);
    // We now know how many processors each unique key maps to. Use this
    // information to create a complete a  new values list, a new keys list, and
    // a processor location list
    std::multimap<int,int> keyMap;
    for (i=0; i<base_keys.size(); i++) {
      keyMap.insert(std::pair<int,int>(base_keys[i],procLoc[i]));
    }
    std::vector<_bus_data_type*> newValues;
    std::vector<int> newKeys;
    std::vector<int> destProcs;
    std::multimap<int,int>::iterator itk;
    _bus_data_type *dptr;
    for (i=0; i<values.size(); i++) {
      itk = keyMap.find(keys[i]);
      if (itk != keyMap.end()) {
        while (itk != keyMap.upper_bound(keys[i])) {
          dptr = new _bus_data_type[nvals];
          for (j=0; j<nvals; j++) {
            dptr[j] = (values[i])[j];
          }
          newValues.push_back(dptr);
          newKeys.push_back(keys[i]);
          destProcs.push_back(itk->second);
          itk++;
        }
      }
    }

    // sort keys into linked list based on which processors they are going to
    int nsize = newKeys.size();
    int ltop[nprocs];
    int ldest[nsize];
    int destNum[nprocs];
    for (i=0; i<nprocs; i++) {
      ltop[i] = -1;
      destNum[i] = 0;
    }
    for (i=0; i<nsize; i++) {
      ldest[i] = -1;
    }
    int iproc;
    for (i=0; i<nsize; i++) {
      iproc = destProcs[i];
      destNum[iproc]++;
      ldest[i] = ltop[iproc]; 
      ltop[iproc] = i;
    }
#ifdef HASH_WITH_MPI
    // Use all-to-all call to determine number of unique keys coming from each
    // processor
    int srcNum[nprocs];
    int ierr;
    int k;
    int one = 1;
    MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
    ierr = MPI_Alltoall(destNum,one,MPI_INT,srcNum,one,MPI_INT,comm);

    // Each process now knows how much data it will receive. Pack data data into
    // an appropriate sized buffer and send it to processors using a all-to-all
    // vector call
    char *sendBuf;
    sendBuf = new char[newValues.size()
      * (sizeof(_bus_data_type)*nvals+sizeof(int))];
    // Fill sendBuf with data from newValues and newKeys
    char *ptr = sendBuf;
    for (i=0; i<nprocs; i++) {
      j = ltop[i];
      if (j>=0) {
        while(j >= 0) {
          ((int*)ptr)[0] = newKeys[j];
          ptr += sizeof(int);
          dptr = (_bus_data_type*)ptr;
          for (k=0; k<nvals; k++) {
            dptr[k] = (newValues[j])[k];
          }
          j = ldest[j];
          ptr += nvals*sizeof(_bus_data_type);
        }
      }
    }

    // Evaluate offsets for data being sent and received and then correct all
    // sizes by the number of bytes in bus_data_pair structure
    int elemsize =  nvals*sizeof(_bus_data_type)+sizeof(int);
    int destOffset[nprocs];
    int srcOffset[nprocs];
    destOffset[0] = 0;
    srcOffset[0] = 0;
    for (i=1; i<nprocs; i++) {
      destOffset[i] = destOffset[i-1] + destNum[i-1];
      srcOffset[i] = srcOffset[i-1] + srcNum[i-1];
    }
    int nvalues = 0;
    for (i=0; i<nprocs; i++) {
      nvalues += srcNum[i];
      destOffset[i] = elemsize*destOffset[i];
      destNum[i] = elemsize*destNum[i];
      srcOffset[i] = elemsize*srcOffset[i];
      srcNum[i] = elemsize*srcNum[i];
    }

    // Allocate buffer to receive data from other processors
    char *recvBuf;
    recvBuf = new char[nvalues*(sizeof(_bus_data_type)*nvals+sizeof(int))];
    
    // Transmit data and clean up buffers that are no longer needed
    ierr = MPI_Alltoallv(sendBuf, destNum, destOffset, MPI_BYTE, recvBuf,
        srcNum, srcOffset, MPI_BYTE, comm);
    delete [] sendBuf;

    // Data is now available on processor that can use it. Pack it into final
    // data structures. Start by creating a map between global and local bus
    // indices on this processors
    keys.clear();
    int vsize = values.size();
    for (i=0; i<vsize; i++) {
      delete [] values[i];
    }
    values.clear();
    int nbus = p_network->numBuses();
    std::multimap<int,int> idxMap;
    std::pair<int,int> idxPair;
    // Create map between global index and local index for locally held buses
    for (i=0; i<nbus; i++) {
      idxPair = std::pair<int,int>(p_network->getOriginalBusIndex(i),i);
      idxMap.insert(idxPair);
    }
    // Get data from recvBuf and pack keys and values arrays
    // with local indices and data
    std::multimap<int,int>::iterator it;
    _bus_data_type *rptr;
    int idx;
    ptr = recvBuf;
    for (i=0; i<nvalues; i++) {
      idx = ((int*)ptr)[0];
      ptr += sizeof(int);
      it = idxMap.find(idx);
      rptr = (_bus_data_type*)ptr;
      ptr += nvals*sizeof(_bus_data_type);
      if (it != idxMap.end()) {
        while (it != idxMap.upper_bound(idx)) {
          dptr = new _bus_data_type[nvals];
          for (j=0; j<nvals; j++) {
            dptr[j] = rptr[j];
          }
          keys.push_back(it->second);
          values.push_back(dptr);
          it++;
        }
      } else {
        printf("p[%d] Unresolved original bus index: %d\n",me,
            idx);
      }
    }

    delete [] recvBuf;
#else
    // Find total number of values going to each processor
    // and get offsets on each processor for each set  of values
    int g_numValues = GA_Create_handle();
    int dims;
    dims = nprocs;
    int one = 1;
    int blocks;
    blocks = 1;
    GA_Set_data(g_numValues,one,&dims,C_INT);
    GA_Set_chunk(g_numValues, &blocks);
    GA_Set_pgroup(g_numValues, p_GAgrp);
    GA_Allocate(g_numValues);
    GA_Zero(g_numValues);
    int r_offset[nprocs];
    for (j=0; j<nprocs; j++) {
      i = (j+me)%nprocs;
      if (destNum[i] > 0) {
        r_offset[i] = NGA_Read_inc(g_numValues,&i,destNum[i]);
      } else {
        r_offset[i] = 0;
      }
    }
    GA_Pgroup_sync(p_GAgrp);
    // Copy values from global array to local array
    int numValues[nprocs];
    int lo, hi;
    lo = 0;
    hi = nprocs-1;
    if (me == 0) {
      if (lo<=hi) NGA_Get(g_numValues,&lo,&hi,numValues,&one);
    } else {
      for (i=0; i<nprocs; i++) {
        numValues[i] = 0;
      }
    }
    char plus[2];
    strcpy(plus,"+");
    GA_Pgroup_igop(p_GAgrp,numValues,nprocs,plus);
    GA_Destroy(g_numValues);

    // Create a global array that can hold all values. Partition the array so
    // that each processor has enough local memory to hold all values coming to
    // it
    int dataSize = nvals*sizeof(_bus_data_type)+sizeof(int);
    int dtype = NGA_Register_type(dataSize);
    // total number of values on all processors
    int totalVals = 0;
    for (i=0; i<nprocs; i++) {
      r_offset[i] += totalVals;
      totalVals += numValues[i];
    }
    if (totalVals == 0) {
      NGA_Deregister_type(dtype);
      return;
    }
    int g_data = GA_Create_handle();
    dims = totalVals;
    GA_Set_data(g_data, one, &dims, dtype);
    GA_Set_pgroup(g_data,p_GAgrp);
    int mapc[nprocs];
    mapc[0] = 0;
    for (i=1; i<nprocs; i++) {
      mapc[i] = mapc[i-1] + numValues[i-1];
    }
    blocks = nprocs;
    GA_Set_irreg_distr(g_data,mapc,&blocks);
    GA_Allocate(g_data);

    // Repack values and send them to the processor with the corresponding buses
    char *bus_data, *ptr;
    int k;
    for (i=0; i<nprocs; i++) {
      j = ltop[i];
      if (j >= 0) {
        bus_data = new char[dataSize*destNum[i]];
        ptr = bus_data;
        while (j >= 0) {
          ((int*)ptr)[0] = newKeys[j];
          ptr += sizeof(int);
          dptr = (_bus_data_type*)ptr;
          for (k=0; k<nvals; k++) {
             dptr[k] = (newValues[j])[k];
          }
          ptr += nvals*sizeof(_bus_data_type);
          j = ldest[j];
        }
        lo = r_offset[i];
        hi = lo + destNum[i] - 1;
        if (lo<=hi) NGA_Put(g_data,&lo,&hi,bus_data,&one);
        delete [] bus_data;
      }
    }

    GA_Pgroup_sync(p_GAgrp);
    // Data is now on processor. Unpack it and put it in arrays for export
    keys.clear();
    int vsize = values.size();
    for (i=0; i<vsize; i++) {
      delete [] values[i];
    }
    values.clear();
    int nbus = p_network->numBuses();
    std::multimap<int,int> idxMap;
    std::pair<int,int> idxPair;
    // Create map between global index and local index for locally held buses
    for (i=0; i<nbus; i++) {
      idxPair = std::pair<int,int>(p_network->getOriginalBusIndex(i),i);
      idxMap.insert(idxPair);
    }
    // Get data from global array and pack keys and values arrays
    // with local indices and data
    int ndata = numValues[me];
    lo = mapc[me];
    hi = lo + ndata - 1;
    if (lo<=hi) NGA_Access(g_data,&lo,&hi,&bus_data,&one);
    std::multimap<int,int>::iterator it;
    int idx;
    ptr = bus_data;
    _bus_data_type *rptr;
    for (i=0; i<ndata; i++) {
      idx = ((int*)ptr)[0];
      ptr += sizeof(int);
      it = idxMap.find(idx);
      rptr = (_bus_data_type*)ptr;
      ptr += nvals*sizeof(_bus_data_type);
      if (it != idxMap.end()) {
        while (it != idxMap.upper_bound(idx)) {
          keys.push_back(it->second);
          dptr = new _bus_data_type[nvals];
          for (j=0; j<nvals; j++) {
            dptr[j] = rptr[j];
          }
          values.push_back(dptr);
          it++;
        }
      } else {
        printf("p[%d] Unresolved original bus index: %d\n",me,
            idx);
      }
    }
    if (lo<=hi) NGA_Release(g_data,&lo,&hi);
    GA_Destroy(g_data);
#endif
#endif
  }

  // Send values corresponding to keys to the processors that own them. On
  // completion, keys contain the local indices of the branches recieving data
  // and values contains the data
  // @param keys on input, a list of integer pair keys corresponding to the
  // branches that receive the data
  // @param branch_ids on output, contains the local indices of branches that
  // receive data
  // @param values on input, a list of values corresponding to the list of keys,
  // on output, the values of the received data
  void distributeBranchValues(std::vector<std::pair<int,int> > &keys,
      std::vector<int> &branch_ids,
      std::vector<_branch_data_type> &values)
  {
#ifdef SYSTOLIC
    int ksize = keys.size();
    int vsize = values.size();
    int me = GA_Pgroup_nodeid(p_GAgrp);
    int nprocs = GA_Pgroup_nnodes(p_GAgrp);
    if (vsize != ksize) {
      char buf[256];
      sprintf(buf,"p[%d] HashDistribution::distributeBranchValues ERROR: length"
          " of keys and values arrays don't match ksize: %d vsize: %d\n",
          me,ksize,vsize);
      printf("%s",buf);
      throw gridpack::Exception(buf);
    }
    
    // Construct mapc array for irregular distribution
    int i;
    int *sizes = new int[nprocs];
    for (i=0; i<nprocs; i++) {
      sizes[i] = 0;
    }
    sizes[me] = ksize;
    char plus[2];
    strcpy(plus,"+");
    GA_Pgroup_igop(p_GAgrp,sizes,nprocs,plus);
    int *mapc = new int[nprocs];
    mapc[0] = 0;
    int total_values = sizes[0];
    for (i=1; i<nprocs; i++) {
      mapc[i] = mapc[i-1]+sizes[i-1];
      total_values += sizes[i];
    }
    if (total_values == 0) {
      delete [] sizes;
      delete [] mapc;
      return;
    }

    // Copy values to a straight array
    branch_data_pair *list;
    if (ksize > 0) {
      list = new branch_data_pair[ksize];
      for (i=0; i<ksize; i++) {
        list[i].idx1 = keys[i].first;
        list[i].idx2 = keys[i].second;
        list[i].data = values[i];
      }
    }
    int g_type = NGA_Register_type(p_size_branch_data);

    // Create global array and store all values in it
    int lo, hi;
    lo = 0;
    for (i=0; i<me; i++) {
      lo += sizes[i];
    }
    hi = lo + sizes[me] - 1;
    int one = 1;
    int g_vals = GA_Create_handle();
    GA_Set_data(g_vals,one,&total_values,g_type);
//    GA_Set_irreg_distr(g_vals,mapc,&nprocs);
    GA_Set_pgroup(g_vals,p_GAgrp);
    if (!GA_Allocate(g_vals)) {
      char buf[256];
      sprintf(buf,"HashDistribution::distributeBranchValues: Unable to allocate"
          " distributed array for storing values");
      printf("%s",buf);
      throw gridpack::Exception(buf);
    }
    if (lo <= hi) NGA_Put(g_vals, &lo, &hi, list, &one);
    GA_Pgroup_sync(p_GAgrp);
    NGA_Deregister_type(g_type);
    if (ksize > 0) delete [] list;
    delete [] mapc;
    delete [] sizes;

    // Create a map that maps original branch indices to local indices
    int idx,idx1,idx2;
    int nbranch = p_network->numBranches();
    std::multimap<std::pair<int,int>,int> hmap;
    for (i=0; i<nbranch; i++) {
      p_network->getOriginalBranchEndpoints(i,&idx1,&idx2);
      hmap.insert(std::pair<std::pair<int,int>,int>(std::pair<int,int>(idx1,idx2),i));
    }
    std::multimap<std::pair<int,int>,int>::iterator it;
    
    // Get values from global array and copy them into output arrays
    branch_ids.clear();
    values.clear();
    double delta = static_cast<double>(total_values)/static_cast<double>(nprocs);
    for (i = 0; i<nprocs; i++) {
      idx = (i+me)%nprocs;
      lo = static_cast<int>(delta*static_cast<double>(idx));
      if (idx<nprocs-1) {
        hi = static_cast<int>(delta*static_cast<double>(idx+1))-1;
      } else {
        hi = total_values-1;
      }
      int nsize = hi - lo + 1;
      if (lo <= hi) {
        list = new branch_data_pair[nsize];
        if (lo<=hi) NGA_Get(g_vals, &lo, &hi, list, &one);
        int j;
        std::pair<int,int> key;
        for (j=0; j<nsize; j++) {
          key = std::pair<int,int>(list[j].idx1,list[j].idx2);
          it = hmap.find(key);
          if (it != hmap.end()) {
            while(it != hmap.upper_bound(key)) {
              branch_ids.push_back(it->second);
              values.push_back(list[j].data);
              it++;
            }
          }
        }
        delete [] list;
      }
    }
    GA_Destroy(g_vals);
#else
    int nprocs = p_network->communicator().size();
    int me = p_network->communicator().rank();
    // Get global indices corresponding to keys
    int i,j;
    // Create a list of unique keys
    std::vector<std::pair<int,int> > base_keys;
    std::set<std::pair<int,int> > key_check;
    std::set<std::pair<int,int> >::iterator itc;
    for (i=0; i<keys.size(); i++) {
      itc = key_check.find(keys[i]);
      if (itc == key_check.end()) {
        key_check.insert(keys[i]);
        base_keys.push_back(keys[i]);
      }
    }
    std::vector<int> procLoc;
    // Get processor that owns global index
    p_indexHashMap->getValues(base_keys,procLoc);
    // We now know how many processors each unique key maps to. Use this
    // information to create a complete a  new values list, a new keys list, and
    // a processor location list
    std::multimap<std::pair<int,int>,int> keyMap;
    for (i=0; i<base_keys.size(); i++) {
      keyMap.insert(std::pair<std::pair<int,int>,int>(base_keys[i],procLoc[i]));
    }
    std::vector<_branch_data_type> newValues;
    std::vector<std::pair<int,int> > newKeys;
    std::vector<int> destProcs;
    j = 0;
    std::multimap<std::pair<int,int>,int>::iterator itk;
    for (i=0; i<values.size(); i++) {
      itk = keyMap.find(keys[i]);
      if (itk != keyMap.end()) {
        while (itk != keyMap.upper_bound(keys[i])) {
          newValues.push_back(values[i]);
          newKeys.push_back(keys[i]);
          destProcs.push_back(itk->second);
          itk++;
        }
      }
    }

    // sort keys into linked list based on which processors they are going to
    int nsize = newKeys.size();
    int ltop[nprocs];
    int ldest[nsize];
    int destNum[nprocs];
    for (i=0; i<nprocs; i++) {
      ltop[i] = -1;
      destNum[i] = 0;
    }
    for (i=0; i<nsize; i++) {
      ldest[i] = -1;
    }
    int iproc;
    for (i=0; i<nsize; i++) {
      iproc = destProcs[i];
      destNum[iproc]++;
      ldest[i] = ltop[iproc]; 
      ltop[iproc] = i;
    }
#ifdef HASH_WITH_MPI
    // Use all-to-all call to determine number of unique keys coming from each
    // processor
    int srcNum[nprocs];
    int ierr;
    int one = 1;
    MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
    ierr = MPI_Alltoall(destNum,one,MPI_INT,srcNum,one,MPI_INT,comm);

    // Each process now knows how much data it will receive. Pack data data into
    // an appropriate sized buffer and send it to processors using a all-to-all
    // vector call
    branch_data_pair *sendBuf;
    sendBuf = new branch_data_pair[newValues.size()];
    // Fill sendBuf with data from newValues and newKeys
    int icnt = 0;
    for (i=0; i<nprocs; i++) {
      j = ltop[i];
      if (j>=0) {
        while(j >= 0) {
          sendBuf[icnt].idx1 = newKeys[j].first;
          sendBuf[icnt].idx2 = newKeys[j].second;
          sendBuf[icnt].data = newValues[j];
          j = ldest[j];
          icnt++;
        }
      }
    }

    // Evaluate offsets for data being sent and received and then correct all
    // sizes by the number of bytes in branch_data_pair structure
    int elemsize =  sizeof(branch_data_pair);
    int destOffset[nprocs];
    int srcOffset[nprocs];
    destOffset[0] = 0;
    srcOffset[0] = 0;
    for (i=1; i<nprocs; i++) {
      destOffset[i] = destOffset[i-1] + destNum[i-1];
      srcOffset[i] = srcOffset[i-1] + srcNum[i-1];
    }
    int nvalues = 0;
    for (i=0; i<nprocs; i++) {
      nvalues += srcNum[i];
      destOffset[i] = elemsize*destOffset[i];
      destNum[i] = elemsize*destNum[i];
      srcOffset[i] = elemsize*srcOffset[i];
      srcNum[i] = elemsize*srcNum[i];
    }

    // Allocate buffer to receive data from other processors
    branch_data_pair *recvBuf;
    recvBuf = new branch_data_pair[nvalues];
    
    // Transmit data and clean up buffers that are no longer needed
    ierr = MPI_Alltoallv(sendBuf, destNum, destOffset, MPI_BYTE, recvBuf,
        srcNum, srcOffset, MPI_BYTE, comm);
    delete [] sendBuf;

    // Data is now available on processor that can use it. Pack it into final
    // data structures. Start by creating a map between global and local branch
    // indices on this processor
    values.clear();
    branch_ids.clear();
    int nbranch = p_network->numBranches();
    std::multimap<std::pair<int,int>,int> idxMap;
    std::pair<std::pair<int,int>,int> idxPair;
    // Create map between global index and local index for locally held buses
    int idx1, idx2;
    std::pair<int,int> key;
    for (i=0; i<nbranch; i++) {
      p_network->getOriginalBranchEndpoints(i,&idx1,&idx2);
      key = std::pair<int,int>(idx1,idx2);
      idxPair = std::pair<std::pair<int,int>,int>(key,i);
      idxMap.insert(idxPair);
    }
    // Get data from recvBuf and pack keys and values arrays
    // with local indices and data
    std::multimap<std::pair<int,int>,int>::iterator it;
    for (i=0; i<nvalues; i++) {
      key = std::pair<int,int>(recvBuf[i].idx1,recvBuf[i].idx2);
      it = idxMap.find(key);
      if (it != idxMap.end()) {
        while (it != idxMap.upper_bound(key)) {
          branch_ids.push_back(it->second);
          values.push_back(recvBuf[i].data);
          it++;
        }
      } else {
        printf("p[%d] Unresolved original branch index: < %d, %d>\n",me,
            recvBuf[i].idx1,recvBuf[i].idx2);
      }
    }

    delete [] recvBuf;
#else
    // Find total number of values going to each processor
    // and get offsets on each processor for each set  of values
    int g_numValues = GA_Create_handle();
    int dims;
    dims = nprocs;
    int one = 1;
    int blocks;
    blocks = 1;
    GA_Set_data(g_numValues,one,&dims,C_INT);
    GA_Set_chunk(g_numValues, &blocks);
    GA_Set_pgroup(g_numValues, p_GAgrp);
    GA_Allocate(g_numValues);
    GA_Zero(g_numValues);
    int r_offset[nprocs];
    for (j=0; j<nprocs; j++) {
      i = (j+me)%nprocs;
      if (destNum[i] > 0) {
        r_offset[i] = NGA_Read_inc(g_numValues,&i,destNum[i]);
      } else {
        r_offset[i] = 0;
      }
    }
    GA_Pgroup_sync(p_GAgrp);
    // Copy values from global array to local array
    int numValues[nprocs];
    int lo, hi;
    lo = 0;
    hi = nprocs-1;
    if (me == 0) {
      if (lo<=hi) NGA_Get(g_numValues,&lo,&hi,numValues,&one);
    } else {
      for (i=0; i<nprocs; i++) {
        numValues[i] = 0;
      }
    }
    char plus[2];
    strcpy(plus,"+");
    GA_Pgroup_igop(p_GAgrp,numValues,nprocs,plus);
    GA_Destroy(g_numValues);

    // Create a global array that can hold all values. Partition the array so
    // that each processor has enough local memory to hold all values coming to
    // it
    int dtype = NGA_Register_type(sizeof(branch_data_pair));
    // total number of values on all processors
    int totalVals = 0;
    for (i=0; i<nprocs; i++) {
      r_offset[i] += totalVals;
      totalVals += numValues[i];
    }
    if (totalVals == 0) {
      NGA_Deregister_type(dtype);
      return;
    }
    int g_data = GA_Create_handle();
    dims = totalVals;
    GA_Set_data(g_data, one, &dims, dtype);
    GA_Set_pgroup(g_data,p_GAgrp);
    int mapc[nprocs];
    mapc[0] = 0;
    for (i=1; i<nprocs; i++) {
      mapc[i] = mapc[i-1] + numValues[i-1];
    }
    blocks = nprocs;
    GA_Set_irreg_distr(g_data,mapc,&blocks);
    GA_Allocate(g_data);

    // Repack values and send them to the processor with the corresponding
    // branches
    branch_data_pair *branch_data;
    int ncnt;
    for (i=0; i<nprocs; i++) {
      j = ltop[i];
      ncnt = 0;
      if (j >= 0) {
        branch_data = new branch_data_pair[destNum[i]];
        while (j >= 0) {
          branch_data[ncnt].idx1 = newKeys[j].first;
          branch_data[ncnt].idx2 = newKeys[j].second;
          branch_data[ncnt].data = newValues[j];
          j = ldest[j];
          ncnt++;
        }
        lo = r_offset[i];
        hi = lo + destNum[i] - 1;
        if (lo<=hi) NGA_Put(g_data,&lo,&hi,branch_data,&one);
        delete [] branch_data;
      }
    }

    GA_Pgroup_sync(p_GAgrp);
    // Data is now on processor. Unpack it and put it in arrays for export
    keys.clear();
    values.clear();
    int nbranch = p_network->numBranches();
    std::multimap<std::pair<int,int>,int> idxMap;
    std::pair<std::pair<int,int>,int> idxPair;
    // Create map between global index and local index for locally held branches
    int idx1, idx2;
    for (i=0; i<nbranch; i++) {
      p_network->getOriginalBranchEndpoints(i,&idx1,&idx2);
      idxPair = std::pair<std::pair<int,int>,int>(std::pair<int,int>(idx1,idx2),i);
      idxMap.insert(idxPair);
    }
    // Get data from global array and pack keys and values arrays
    // with local indices and data
    int ndata = numValues[me];
    lo = mapc[me];
    hi = lo + ndata - 1;
    if (lo<=hi) NGA_Access(g_data,&lo,&hi,&branch_data,&one);
    std::multimap<std::pair<int,int>,int>::iterator it;
    std::pair<int,int> key;
    for (i=0; i<ndata; i++) {
      key = std::pair<int,int>(branch_data[i].idx1,branch_data[i].idx2);
      it = idxMap.find(key);
      if (it != idxMap.end()) {
        while (it != idxMap.upper_bound(key)) {
          branch_ids.push_back(it->second);
          values.push_back(branch_data[i].data);
          it++;
        }
      } else {
        printf("p[%d] Unresolved original branch index: < %d, %d >\n",me,
            branch_data[i].idx1,branch_data[i].idx2);
      }
    }
    if (lo<=hi) NGA_Release(g_data,&lo,&hi);
    GA_Destroy(g_data);
#endif
#endif
  }

  // Send values corresponding to keys to the processors that own them. On
  // completion, branch_ids contain the local indices of the branches recieving data
  // and values contains the data
  // @param keys on input, a list of integer pair keys corresponding to the
  // branches that receive the data
  // @param branch_ids on output, contains local indices of branches that
  // receive data
  // @param values on input, a list of values corresponding to the list of keys,
  // on output, the values of the received data
  // @param nvals the number of values in each array
  void distributeBranchValues(std::vector<std::pair<int,int> > &keys,
      std::vector<int> &branch_ids,
      std::vector<_branch_data_type*> &values, int nvals)
  {
#ifdef SYSTOLIC
    int ksize = keys.size();
    int vsize = values.size();
    int me = GA_Pgroup_nodeid(p_GAgrp);
    int nprocs = GA_Pgroup_nnodes(p_GAgrp);
    if (vsize != ksize) {
      char buf[256];
      sprintf(buf,"p[%d] HashDistribution::distributeBranchValues ERROR: length"
          " of keys and values arrays don't match ksize: %d vsize: %d\n",
          me,ksize,vsize);
      printf("%s",buf);
      throw gridpack::Exception(buf);
    }
    
    // Construct mapc array for irregular distribution
    int i, j;
    int *sizes = new int[nprocs];
    for (i=0; i<nprocs; i++) {
      sizes[i] = 0;
    }
    sizes[me] = ksize;
    char plus[2];
    strcpy(plus,"+");
    GA_Pgroup_igop(p_GAgrp,sizes,nprocs,plus);
    int *mapc = new int[nprocs];
    mapc[0] = 0;
    int total_values = sizes[0];
    for (i=1; i<nprocs; i++) {
      mapc[i] = mapc[i-1]+sizes[i-1];
      total_values += sizes[i];
    }
    if (total_values == 0) {
      delete [] sizes;
      delete [] mapc;
      return;
    }

    // Copy values to a straight array
    char *list, *ptr;
    if (ksize > 0) {
      list = new char[ksize*(nvals*sizeof(_branch_data_type)+2*sizeof(int))];
      ptr = list;
      for (i=0; i<ksize; i++) {
        ((int*)ptr)[0] = keys[i].first;
        ((int*)ptr)[1] = keys[i].second;
        ptr += 2*sizeof(int);
        for (j=0; j<nvals; j++) {
          ((_branch_data_type*)ptr)[j] = (values[i])[j];
        }
        ptr += nvals*sizeof(_branch_data_type);
      }
    }
    int g_type = NGA_Register_type(nvals*sizeof(_branch_data_type)+2*sizeof(int));

    // Create global array and store all values in it
    int lo, hi;
    lo = 0;
    for (i=0; i<me; i++) {
      lo += sizes[i];
    }
    hi = lo + sizes[me] - 1;
    int one = 1;
    int g_vals = GA_Create_handle();
    GA_Set_data(g_vals,one,&total_values,g_type);
//    GA_Set_irreg_distr(g_vals,mapc,&nprocs);
    GA_Set_pgroup(g_vals,p_GAgrp);
    if (!GA_Allocate(g_vals)) {
      char buf[256];
      sprintf(buf,"HashDistribution::distributeBranchValues: Unable to allocate"
          " distributed array for storing values");
      printf("%s",buf);
      throw gridpack::Exception(buf);
    }
    if (lo <= hi) NGA_Put(g_vals, &lo, &hi, list, &one);
    GA_Pgroup_sync(p_GAgrp);
    NGA_Deregister_type(g_type);
    if (ksize > 0) delete [] list;
    delete [] mapc;
    delete [] sizes;

    // Create a map that maps original branch indices to local indices
    int idx,idx1,idx2;
    int nbranch = p_network->numBranches();
    std::multimap<std::pair<int,int>,int> hmap;
    for (i=0; i<nbranch; i++) {
      p_network->getOriginalBranchEndpoints(i,&idx1,&idx2);
      hmap.insert(std::pair<std::pair<int,int>,int>(std::pair<int,int>(idx1,idx2),i));
    }
    std::multimap<std::pair<int,int>,int>::iterator it;
    
    // Get values from global array and copy them into output arrays
    branch_ids.clear();
    for (i=0; i<ksize; i++) {
      delete [] values[i];
    }
    values.clear();
    double delta = static_cast<double>(total_values)/static_cast<double>(nprocs);
    for (i = 0; i<nprocs; i++) {
      idx = (i+me)%nprocs;
      lo = static_cast<int>(delta*static_cast<double>(idx));
      if (idx<nprocs-1) {
        hi = static_cast<int>(delta*static_cast<double>(idx+1))-1;
      } else {
        hi = total_values-1;
      }
      int nsize = hi - lo + 1;
      if (lo <= hi) {
        list = new char[nsize*(nvals*sizeof(_branch_data_type)+2*sizeof(int))];
        if (lo<=hi) NGA_Get(g_vals, &lo, &hi, list, &one);
        int j, k;
        std::pair<int,int> key;
        ptr = list;
        for (j=0; j<nsize; j++) {
          key = std::pair<int,int>(((int*)ptr)[0],((int*)ptr)[1]);
          ptr += 2*sizeof(int);
          it = hmap.find(key);
          if (it != hmap.end()) {
            while(it != hmap.upper_bound(key)) {
              branch_ids.push_back(it->second);
              _branch_data_type *data = new _branch_data_type[nvals];
              for (k=0; k<nvals; k++) {
                data[k] = ((_branch_data_type*)ptr)[k];
              }
              values.push_back(data);
              it++;
            }
          }
          ptr += nvals*sizeof(_branch_data_type);
        }
        delete [] list;
      }
    }
    GA_Destroy(g_vals);
#else
    int nprocs = p_network->communicator().size();
    int me = p_network->communicator().rank();
    // Get global indices corresponding to keys
    int i,j;
    // Create a list of unique keys
    std::vector<std::pair<int,int> > base_keys;
    std::set<std::pair<int,int> > key_check;
    std::set<std::pair<int,int> >::iterator itc;
    for (i=0; i<keys.size(); i++) {
      itc = key_check.find(keys[i]);
      if (itc == key_check.end()) {
        key_check.insert(keys[i]);
        base_keys.push_back(keys[i]);
      }
    }
    std::vector<int> procLoc;
    // Get processor that owns global index
    p_indexHashMap->getValues(base_keys,procLoc);
    // We now know how many processors each unique key maps to. Use this
    // information to create a complete a  new values list, a new keys list, and
    // a processor location list
    std::multimap<std::pair<int,int>,int> keyMap;
    for (i=0; i<base_keys.size(); i++) {
      keyMap.insert(std::pair<std::pair<int,int>,int>(base_keys[i],procLoc[i]));
    }
    std::vector<_branch_data_type*> newValues;
    std::vector<std::pair<int,int> > newKeys;
    std::vector<int> destProcs;
    std::multimap<std::pair<int,int>,int>::iterator itk;
    _branch_data_type *dptr;
    for (i=0; i<values.size(); i++) {
      itk = keyMap.find(keys[i]);
      if (itk != keyMap.end()) {
        while (itk != keyMap.upper_bound(keys[i])) {
          dptr = new _branch_data_type[nvals];
          for (j=0; j<nvals; j++) {
            dptr[j] = (values[i])[j];
          }
          newValues.push_back(values[i]);
          newKeys.push_back(keys[i]);
          destProcs.push_back(itk->second);
          itk++;
        }
      }
    }

    // sort keys into linked list based on which processors they are going to
    int nsize = newKeys.size();
    int ltop[nprocs];
    int ldest[nsize];
    int destNum[nprocs];
    for (i=0; i<nprocs; i++) {
      ltop[i] = -1;
      destNum[i] = 0;
    }
    for (i=0; i<nsize; i++) {
      ldest[i] = -1;
    }
    int iproc;
    for (i=0; i<nsize; i++) {
      iproc = destProcs[i];
      destNum[iproc]++;
      ldest[i] = ltop[iproc]; 
      ltop[iproc] = i;
    }
#ifdef HASH_WITH_MPI
    // Use all-to-all call to determine number of unique keys coming from each
    // processor
    int srcNum[nprocs];
    int ierr;
    int one = 1;
    MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
    ierr = MPI_Alltoall(destNum,one,MPI_INT,srcNum,one,MPI_INT,comm);

    // Each process now knows how much data it will receive. Pack data data into
    // an appropriate sized buffer and send it to processors using a all-to-all
    // vector call
    char *sendBuf;
    sendBuf = new char[newValues.size()
      * (sizeof(_branch_data_type)*nvals+2*sizeof(int))];
    // Fill sendBuf with data from newValues and newKeys
    char *ptr = sendBuf;
    int k;
    for (i=0; i<nprocs; i++) {
      j = ltop[i];
      if (j>=0) {
        while(j >= 0) {
          ((int*)ptr)[0] = newKeys[j].first;
          ((int*)ptr)[1] = newKeys[j].second;
          ptr += 2*sizeof(int);
          dptr = (_branch_data_type*)ptr;
          for (k=0; k<nvals; k++) {
            dptr[k] = (newValues[j])[k];
          }
          j = ldest[j];
          ptr += nvals*sizeof(_branch_data_type);
        }
      }
    }

    // Evaluate offsets for data being sent and received and then correct all
    // sizes by the number of bytes in branch_data_pair structure
    int elemsize =  nvals*sizeof(_branch_data_type)+2*sizeof(int);
    int destOffset[nprocs];
    int srcOffset[nprocs];
    destOffset[0] = 0;
    srcOffset[0] = 0;
    for (i=1; i<nprocs; i++) {
      destOffset[i] = destOffset[i-1] + destNum[i-1];
      srcOffset[i] = srcOffset[i-1] + srcNum[i-1];
    }
    int nvalues = 0;
    for (i=0; i<nprocs; i++) {
      nvalues += srcNum[i];
      destOffset[i] = elemsize*destOffset[i];
      destNum[i] = elemsize*destNum[i];
      srcOffset[i] = elemsize*srcOffset[i];
      srcNum[i] = elemsize*srcNum[i];
    }

    // Allocate buffer to receive data from other processors
    char *recvBuf;
    recvBuf = new char[nvalues*(sizeof(_branch_data_type)*nvals+2*sizeof(int))];
    
    // Transmit data and clean up buffers that are no longer needed
    ierr = MPI_Alltoallv(sendBuf, destNum, destOffset, MPI_BYTE, recvBuf,
        srcNum, srcOffset, MPI_BYTE, comm);
    delete [] sendBuf;

    // Data is now available on processor that can use it. Pack it into final
    // data structures. Start by creating a map between global and local branch
    // indices on this processor
    int vsize = values.size();
    for (i=0; i<vsize; i++) {
      delete [] values[i];
    }
    values.clear();
    branch_ids.clear();
    int nbranch = p_network->numBranches();
    std::multimap<std::pair<int,int>,int> idxMap;
    std::pair<std::pair<int,int>,int> idxPair;
    // Create map between global index and local index for locally held buses
    int idx1, idx2;
    std::pair<int,int> key;
    for (i=0; i<nbranch; i++) {
      p_network->getOriginalBranchEndpoints(i,&idx1,&idx2);
      key = std::pair<int,int>(idx1,idx2);
      idxPair = std::pair<std::pair<int,int>,int>(key,i);
      idxMap.insert(idxPair);
    }
    // Get data from recvBuf and pack keys and values arrays
    // with local indices and data
    std::multimap<std::pair<int,int>,int>::iterator it;
    _branch_data_type *rptr;
    ptr = recvBuf;
    for (i=0; i<nvalues; i++) {
      key = std::pair<int,int>(((int*)ptr)[0], ((int*)ptr)[1]);
      ptr += 2*sizeof(int);
      it = idxMap.find(key);
      rptr = (_branch_data_type*)ptr;
      ptr += nvals*sizeof(_branch_data_type);
      if (it != idxMap.end()) {
        while (it != idxMap.upper_bound(key)) {
          dptr = new _branch_data_type[nvals];
          for (j=0; j<nvals; j++) {
            dptr[j] = rptr[j];
          }
          branch_ids.push_back(it->second);
          values.push_back(dptr);
          it++;
        }
      } else {
        printf("p[%d] Unresolved original branch index: < %d, %d>\n",me,
            key.first,key.second);
      }
    }

    delete [] recvBuf;
#else
    // Find total number of values going to each processor
    // and get offsets on each processor for each set  of values
    int g_numValues = GA_Create_handle();
    int dims;
    dims = nprocs;
    int one = 1;
    int blocks;
    blocks = 1;
    GA_Set_data(g_numValues,one,&dims,C_INT);
    GA_Set_chunk(g_numValues, &blocks);
    GA_Set_pgroup(g_numValues, p_GAgrp);
    GA_Allocate(g_numValues);
    GA_Zero(g_numValues);
    int r_offset[nprocs];
    for (j=0; j<nprocs; j++) {
      i = (j+me)%nprocs;
      if (destNum[i] > 0) {
        r_offset[i] = NGA_Read_inc(g_numValues,&i,destNum[i]);
      } else {
        r_offset[i] = 0;
      }
    }
    GA_Pgroup_sync(p_GAgrp);
    // Copy values from global array to local array
    int numValues[nprocs];
    int lo, hi;
    lo = 0;
    hi = nprocs-1;
    if (me == 0) {
      if (lo<=hi) NGA_Get(g_numValues,&lo,&hi,numValues,&one);
    } else {
      for (i=0; i<nprocs; i++) {
        numValues[i] = 0;
      }
    }
    char plus[2];
    strcpy(plus,"+");
    GA_Pgroup_igop(p_GAgrp,numValues,nprocs,plus);
    GA_Destroy(g_numValues);

    // Create a global array that can hold all values. Partition the array so
    // that each processor has enough local memory to hold all values coming to
    // it
    int dataSize = nvals*sizeof(_branch_data_type)+2*sizeof(int);
    int dtype = NGA_Register_type(dataSize);
    // total number of values on all processors
    int totalVals = 0;
    for (i=0; i<nprocs; i++) {
      r_offset[i] += totalVals;
      totalVals += numValues[i];
    }
    if (totalVals == 0) {
      NGA_Deregister_type(dtype);
      return;
    }
    int g_data = GA_Create_handle();
    dims = totalVals;
    GA_Set_data(g_data, one, &dims, dtype);
    GA_Set_pgroup(g_data,p_GAgrp);
    int mapc[nprocs];
    mapc[0] = 0;
    for (i=1; i<nprocs; i++) {
      mapc[i] = mapc[i-1] + numValues[i-1];
    }
    blocks = nprocs;
    GA_Set_irreg_distr(g_data,mapc,&blocks);
    GA_Allocate(g_data);

    // Repack values and send them to the processor with the corresponding
    // branches
    char *branch_data, *ptr;
    int k;
    for (i=0; i<nprocs; i++) {
      j = ltop[i];
      if (j >= 0) {
        branch_data = new char[dataSize*destNum[i]];
        ptr = branch_data;
        while (j >= 0) {
          ((int*)ptr)[0] = newKeys[j].first;
          ((int*)ptr)[1] = newKeys[j].second;
          ptr += 2*sizeof(int);
          dptr = (_branch_data_type*)ptr;
          for (k=0; k<nvals; k++) {
            dptr[k] = (newValues[j])[k];
          }
          ptr += nvals*sizeof(_branch_data_type);
          j = ldest[j];
        }
        lo = r_offset[i];
        hi = lo + destNum[i] - 1;
        if (lo<=hi) NGA_Put(g_data,&lo,&hi,branch_data,&one);
        delete [] branch_data;
      }
    }

    GA_Pgroup_sync(p_GAgrp);
    // Data is now on processor. Unpack it and put it in arrays for export
    keys.clear();
    int vsize = values.size();
    for (i=0; i<vsize; i++) {
      delete [] values[i];
    }
    values.clear();
    int nbranch = p_network->numBranches();
    std::multimap<std::pair<int,int>,int> idxMap;
    std::pair<std::pair<int,int>,int> idxPair;
    // Create map between global index and local index for locally held branches
    int idx1, idx2;
    for (i=0; i<nbranch; i++) {
      p_network->getOriginalBranchEndpoints(i,&idx1,&idx2);
      idxPair = std::pair<std::pair<int,int>,int>(std::pair<int,int>(idx1,idx2),i);
      idxMap.insert(idxPair);
    }
    // Get data from global array and pack keys and values arrays
    // with local indices and data
    int ndata = numValues[me];
    lo = mapc[me];
    hi = lo + ndata - 1;
    if (lo<=hi) NGA_Access(g_data,&lo,&hi,&branch_data,&one);
    std::multimap<std::pair<int,int>,int>::iterator it;
    std::pair<int,int> key;
    ptr = branch_data;
    _branch_data_type *rptr;
    for (i=0; i<ndata; i++) {
      key = std::pair<int,int>(((int*)ptr)[0], ((int*)ptr)[1]);
      ptr += 2*sizeof(int);
      it = idxMap.find(key);
      rptr = (_branch_data_type*)ptr;
      ptr += nvals*sizeof(_branch_data_type);
      if (it != idxMap.end()) {
        while (it != idxMap.upper_bound(key)) {
          branch_ids.push_back(it->second);
          dptr = new _branch_data_type[nvals];
          for (j=0; j<nvals; j++) {
            dptr[j] = rptr[j];
          }
          values.push_back(dptr);
          it++;
        }
      } else {
        printf("p[%d] Unresolved original branch index: < %d, %d >\n",me,
            key.first,key.second);
      }
    }
    if (lo<=hi) NGA_Release(g_data,&lo,&hi);
    GA_Destroy(g_data);
#endif
#endif
  }

private:

  // processor(s) that own  original bus index or bus index pair
  boost::shared_ptr<gridpack::hash_map::GlobalIndexHashMap> p_indexHashMap;

  NetworkPtr p_network;

  int p_size_bus_data;
  int p_size_branch_data;

  int p_GAgrp;

};


} // namespace gridpack
} // namespace hash_distr

#endif

