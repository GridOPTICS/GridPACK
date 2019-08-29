/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   index_hash.cpp
 * @author Bruce Palmer
 * @date   2014-06-24 09:25:03 d3g293
 * 
 * @brief  
 * This is a utility that is designed to provide a relatively efficient way of
 * mapping between different sets of indexes in a distributed way. Note that all
 * these operations are collective
 * 
 */

// -------------------------------------------------------------

#include <ga.h>
#include "index_hash.hpp"

//#define HASH_WITH_MPI


// -------------------------------------------------------------
//  class GlobalIndexHashMap
// -------------------------------------------------------------

namespace gridpack {
namespace hash_map {

// Default constructor
GlobalIndexHashMap::GlobalIndexHashMap(const parallel::Communicator &comm)
{
  p_nprocs = comm.size();
  p_me = comm.rank();
  p_comm = static_cast<MPI_Comm>(comm);
  p_GAgrp = comm.getGroup();
}

// Default destructor
GlobalIndexHashMap::~GlobalIndexHashMap(void)
{
}

// add key-value pairs to hash map where key is single integer.
// @param pairs list of key-value pairs where both keys and values are
//              integers
void GlobalIndexHashMap::addPairs(std::vector<std::pair<int,int> > &pairs)
{
  // Need to distribute key-value pairs between processors based on the value
  // returned by the hashValue function. Start by constructing a linked list of
  // where each pair needs to go
  int i, j;
  int size = pairs.size();
  // ndest[idx] number of values to send to processor idx, nrecv[idx] is number
  // of values received from processor idx
  std::vector<int> ndest(p_nprocs);
  std::vector<int> nrecv(p_nprocs);
  // arrays for linked list
  std::vector<int> ltop(p_nprocs);
  std::vector<int> ldest(size);
  // initialize all arrays
  for (i=0; i<p_nprocs; i++) {
    ndest[i] = 0;
    nrecv[i] = 0;
    ltop[i] = -1;
  }
  for (i=0; i<size; i++) {
    ldest[i] = -1;
  }
  // create linked list
  int hash;
  for (i=0; i<size; i++) {
    hash = hashValue(pairs[i].first);
    ndest[hash]++;
    ldest[i] = ltop[hash];
    ltop[hash] = i;
  }
#ifdef HASH_WITH_MPI
  // send data to processors based on linked list. Start by evaluating how much
  // data will be received from other processors using an all-to-all call
  int ierr;
  int one = 1;
  ierr = MPI_Alltoall(&ndest[0], one, MPI_INT, &nrecv[0], one, MPI_INT, p_comm);
  // nrecv now contains the number of pairs that will be received from other
  // processors. Use this information to set up send and receive along with
  // their offsets for a all-to-all-v call
  int rsize = 0;
  for (i=0; i<p_nprocs; i++) {
    rsize += nrecv[i];
    ndest[i] = 2*ndest[i];
    nrecv[i] = 2*nrecv[i];
  }
  std::vector<int> send_pair(2*size);
  std::vector<int> recv_pair(2*rsize);
  std::vector<int> s_offsets(p_nprocs);
  std::vector<int> r_offsets(p_nprocs);
  s_offsets[0] = 0;
  r_offsets[0] = 0;
  for (i=1; i<p_nprocs; i++) {
    s_offsets[i] = s_offsets[i-1]+ndest[i-1];
    r_offsets[i] = r_offsets[i-1]+nrecv[i-1];
  }
  // Fill up send_pair array with key-value pairs using the linked list
  for (i=0; i<p_nprocs; i++) {
    int offset = s_offsets[i];
    if (ndest[i] > 0) {
      int count = offset;
      int j = ltop[i];
      while (j >= 0) {
        send_pair[count] = pairs[j].first;
        send_pair[count+1] = pairs[j].second;
        j = ldest[j];
        count += 2;
      }
    }
  }
  // Send buffer is full so distribute contents to all processors
  ierr = MPI_Alltoallv(&send_pair[0], &ndest[0], &s_offsets[0], MPI_INT,
      &recv_pair[0], &nrecv[0], &r_offsets[0], MPI_INT, p_comm);
  // key-value pairs are available, so set up local hash table
  p_umap.clear();
  int first, second;
  for (i=0; i<rsize; i++) {
    first = recv_pair[2*i];
    second = recv_pair[2*i+1];
    p_umap.insert(std::pair<int, int>(first,second));
  }
#else
  // Create a global array to count how many values are coming from each
  // processor and use this to create a set of offsets
  int g_offset = GA_Create_handle();
  int dims = p_nprocs;
  int one = 1;
  int blocks = 1;
  GA_Set_data(g_offset,one,&dims,C_INT);
  GA_Set_chunk(g_offset, &blocks);
  GA_Set_pgroup(g_offset, p_GAgrp);
  GA_Allocate(g_offset);
  GA_Zero(g_offset);

  // Get offsets on remote processors
  std::vector<int> r_offset(p_nprocs);
  for (j=0; j<p_nprocs; j++) {
    i = (j+p_me)%p_nprocs;
    if (ndest[i] > 0) {
      r_offset[i] = NGA_Read_inc(g_offset, &i, ndest[i]);
    } else {
      r_offset[i] = 0;
    }
  }
  GA_Pgroup_sync(p_GAgrp);
  // Copy values from global array to local array and evaluate global offsets
  std::vector<int> numValues(p_nprocs);
  int lo, hi;
  lo = 0;
  hi = p_nprocs-1;
  if (p_me == 0) {
    NGA_Get(g_offset,&lo,&hi,&numValues[0],&one);
  } else {
    for (i=0; i<p_nprocs; i++) {
      numValues[i] = 0;
    }
  }
  char plus[2];
  strcpy(plus,"+");
  GA_Pgroup_igop(p_GAgrp, &numValues[0], p_nprocs, plus);
  GA_Destroy(g_offset);
  int totalVals = 0;
  for (i=0; i<p_nprocs; i++) {
    r_offset[i] += totalVals;
    totalVals += numValues[i];
  }
  if (totalVals == 0) return;

  // Create a global array that can hold all values
  int dtype = NGA_Register_type(sizeof(std::pair<int,int>));
  int g_data = GA_Create_handle();
  dims = totalVals+1;
  GA_Set_data(g_data,one,&dims,dtype);
  GA_Set_pgroup(g_data,p_GAgrp);
  std::vector<int> mapc(p_nprocs);
  mapc[0] = 0;
  for (i=1; i<p_nprocs; i++) {
    mapc[i] = mapc[i-1] + numValues[i-1];
  }
  blocks = p_nprocs;
  GA_Set_irreg_distr(g_data,&mapc[0],&blocks);
  GA_Allocate(g_data);
  NGA_Deregister_type(dtype);

  // Repack values and send them to the processor that owns the corresponding
  // keys
  std::pair<int,int> *data_pairs;
  int ncnt;
  for (i=0; i<p_nprocs; i++) {
    j = ltop[i];
    ncnt = 0;
    if (j >= 0) {
      data_pairs = new std::pair<int,int>[ndest[i]];
      while (j >= 0) {
        data_pairs[ncnt] = pairs[j];
        j = ldest[j];
        ncnt++;
      }
      lo = r_offset[i];
      hi = lo + ndest[i] - 1;
      NGA_Put(g_data, &lo, &hi, data_pairs, &one);
      delete [] data_pairs;
    }
  }
  GA_Pgroup_sync(p_GAgrp);

  // Data is now on the processor. Store it in hash map
  lo = mapc[p_me];
  hi = mapc[p_me] + numValues[p_me] - 1;
  int ld;
  p_umap.clear();
  if (lo<=hi) {
    NGA_Access(g_data,&lo,&hi,&data_pairs,&ld);
    for (i=0; i<numValues[p_me]; i++) {
      p_umap.insert(data_pairs[i]);
    }
    NGA_Release(g_data,&lo,&hi);
  }
  GA_Destroy(g_data);
#endif
}

// add key-value pairs to hash map where key is another index pair of integers
// @param pairs list of key-value pairs where key is a pair of integers and
//              value is a single integer
void GlobalIndexHashMap::addPairs(std::vector<std::pair<std::pair<int,int>,int> > &pairs)
{
  // Need to distribute key-value pairs between processors based on the value
  // returned by the hashValue function. Start by constructing a linked list of
  // where each pair needs to go
  int i, j;
  int size = pairs.size();
  // ndest[idx] number of values to send to processor idx, nrecv[idx] is number
  // of values received from processor idx
  std::vector<int> ndest(p_nprocs);
  std::vector<int> nrecv(p_nprocs);
  // arrays for linked list
  std::vector<int> ltop(p_nprocs);
  std::vector<int> ldest(size);
  // initialize all arrays
  for (i=0; i<p_nprocs; i++) {
    ndest[i] = 0;
    nrecv[i] = 0;
    ltop[i] = -1;
  }
  for (i=0; i<size; i++) {
    ldest[i] = -1;
  }
  // create linked list
  int hash;
  for (i=0; i<size; i++) {
    hash = pairHashValue(pairs[i].first);
    ndest[hash]++;
    ldest[i] = ltop[hash];
    ltop[hash] = i;
  }
#ifdef HASH_WITH_MPI
  // send data to processors based on linked list. Start by evaluating how much
  // data will be received from other processors using an all-to-all call
  int ierr;
  int one = 1;
  ierr = MPI_Alltoall(&ndest[0], one, MPI_INT, &nrecv[0], one, MPI_INT, p_comm);
  // nrecv now contains the number of pairs that will be received from other
  // processors. Use this information to set up send and receive along with
  // their offsets for a all-to-all-v call
  int rsize = 0;
  for (i=0; i<p_nprocs; i++) {
    rsize += nrecv[i];
    ndest[i] = 3*ndest[i];
    nrecv[i] = 3*nrecv[i];
  }
  std::vector<int> send_pair(3*size);
  std::vector<int> recv_pair(3*rsize);
  std::vector<int> s_offsets(p_nprocs);
  std::vector<int> r_offsets(p_nprocs);
  s_offsets[0] = 0;
  r_offsets[0] = 0;
  for (i=1; i<p_nprocs; i++) {
    s_offsets[i] = s_offsets[i-1]+ndest[i-1];
    r_offsets[i] = r_offsets[i-1]+nrecv[i-1];
  }
  // Fill up send_pair array with key-value pairs using the linked list
  for (i=0; i<p_nprocs; i++) {
    int offset = s_offsets[i];
    if (ndest[i] > 0) {
      int count = offset;
      int j = ltop[i];
      while (j >= 0) {
        send_pair[count] = (pairs[j].first).first;
        send_pair[count+1] = (pairs[j].first).second;
        send_pair[count+2] = pairs[j].second;
        j = ldest[j];
        count += 3;
      }
    }
  }
  // Send buffer is full so distribute contents to all processors
  ierr = MPI_Alltoallv(&send_pair[0], &ndest[0], &s_offsets[0], MPI_INT,
      &recv_pair[0], &nrecv[0], &r_offsets[0], MPI_INT, p_comm);
  // key-value pairs are available, so set up local hash table
  p_pmap.clear();
  int first, second, third;
  for (i=0; i<rsize; i++) {
    first = recv_pair[3*i];
    second = recv_pair[3*i+1];
    third = recv_pair[3*i+2];
    p_pmap.insert(std::pair<std::pair<int,int>,int>(std::pair<int, int>(first,second),third));
  }
#else
  // Create a global array to count how many values are coming from each
  // processor and use this to create a set of offsets
  int g_offset = GA_Create_handle();
  int dims = p_nprocs;
  int one = 1;
  int blocks = 1;
  GA_Set_data(g_offset,one,&dims,C_INT);
  GA_Set_chunk(g_offset, &blocks);
  GA_Set_pgroup(g_offset, p_GAgrp);
  GA_Allocate(g_offset);
  GA_Zero(g_offset);

  // Get offsets on remote processors
  std::vector<int> r_offset(p_nprocs);
  for (j=0; j<p_nprocs; j++) {
    i = (j+p_me)%p_nprocs;
    if (ndest[i] > 0) {
      r_offset[i] = NGA_Read_inc(g_offset, &i, ndest[i]);
    } else {
      r_offset[i] = 0;
    }
  }
  GA_Pgroup_sync(p_GAgrp);
  // Copy values from global array to local array and evaluate global offsets
  std::vector<int> numValues(p_nprocs);
  int lo, hi;
  lo = 0;
  hi = p_nprocs-1;
  if (p_me == 0) {
    NGA_Get(g_offset,&lo,&hi,&numValues[0],&one);
  } else {
    for (i=0; i<p_nprocs; i++) {
      numValues[i] = 0;
    }
  }
  char plus[2];
  strcpy(plus,"+");
  GA_Pgroup_igop(p_GAgrp, &numValues[0], p_nprocs, plus);
  GA_Destroy(g_offset);
  int totalVals = 0;
  for (i=0; i<p_nprocs; i++) {
    r_offset[i] += totalVals;
    totalVals += numValues[i];
  }
  if (totalVals == 0) return;

  // Create a global array that can hold all values
  int dtype = NGA_Register_type(sizeof(std::pair<std::pair<int,int>,int>));
  int g_data = GA_Create_handle();
  dims = totalVals+1;
  GA_Set_data(g_data,one,&dims,dtype);
  GA_Set_pgroup(g_data,p_GAgrp);
  std::vector<int> mapc(p_nprocs);
  mapc[0] = 0;
  for (i=1; i<p_nprocs; i++) {
    mapc[i] = mapc[i-1] + numValues[i-1];
  }
  blocks = p_nprocs;
  GA_Set_irreg_distr(g_data,&mapc[0],&blocks);
  GA_Allocate(g_data);
  NGA_Deregister_type(dtype);

  // Repack values and send them to the processor that owns the corresponding
  // keys
  std::pair<std::pair<int,int>,int> *data_pairs;
  int ncnt;
  for (i=0; i<p_nprocs; i++) {
    j = ltop[i];
    ncnt = 0;
    if (j >= 0) {
      data_pairs = new std::pair<std::pair<int,int>,int>[ndest[i]];
      while (j >= 0) {
        data_pairs[ncnt] = pairs[j];
        j = ldest[j];
        ncnt++;
      }
      lo = r_offset[i];
      hi = lo + ndest[i] - 1;
      NGA_Put(g_data, &lo, &hi, data_pairs, &one);
      delete [] data_pairs;
    }
  }
  GA_Pgroup_sync(p_GAgrp);

  // Data is now on the processor. Store it in hash map
  lo = mapc[p_me];
  hi = mapc[p_me] + numValues[p_me] - 1;
  int ld;
  p_pmap.clear();
  if (lo<=hi) {
    NGA_Access(g_data,&lo,&hi,&data_pairs,&ld);
    for (i=0; i<numValues[p_me]; i++) {
      p_pmap.insert(data_pairs[i]);
    }
    NGA_Release(g_data,&lo,&hi);
  }
  GA_Destroy(g_data);
#endif
}

// get values corresponding to a list of keys from the hash map where key is a
// single integer
// @param keys list of integer keys
// @param values returned list of values corresponding to the list of keys
void GlobalIndexHashMap::getValues(std::vector<int> &keys, std::vector<int> &values)
{
  // Need to distribute keys to processors that hold the corresponding values.
  // Start by constructing a linked list of where each key needs to go
  int i, j;
  int size = keys.size();
  // ndest[idx] number of values to send to processor idx, nrecv[idx] is number
  // of values received from processor idx
  std::vector<int> ndest(p_nprocs);
  std::vector<int> nrecv(p_nprocs);
  // arrays for linked list
  std::vector<int> ltop(p_nprocs);
  std::vector<int> ldest(size);
  // initialize all arrays
  for (i=0; i<p_nprocs; i++) {
    ndest[i] = 0;
    nrecv[i] = 0;
    ltop[i] = -1;
  }
  for (i=0; i<size; i++) {
    ldest[i] = -1;
  }
  // create linked list
  int hash;
  for (i=0; i<size; i++) {
    hash = hashValue(keys[i]);
    ndest[hash]++;
    ldest[i] = ltop[hash];
    ltop[hash] = i;
  }
#ifdef HASH_WITH_MPI
  // send keys to processors based on linked list. Start by evaluating how much
  // data will be received from other processors using an all-to-all call
  int ierr;
  int one = 1;
  ierr = MPI_Alltoall(&ndest[0], one, MPI_INT, &nrecv[0], one, MPI_INT, p_comm);
  // nrecv now contains the number of keys that will be received from other
  // processors. Use this information to set up send and receive along with
  // their offsets for an all-to-all-v call
  int rsize = 0;
  for (i=0; i<p_nprocs; i++) {
    rsize += nrecv[i];
  }
  std::vector<int> send_keys(size);
  std::vector<int> recv_keys(rsize);
  std::vector<int> s_offsets(p_nprocs);
  std::vector<int> r_offsets(p_nprocs);
  s_offsets[0] = 0;
  r_offsets[0] = 0;
  for (i=1; i<p_nprocs; i++) {
    s_offsets[i] = s_offsets[i-1]+ndest[i-1];
    r_offsets[i] = r_offsets[i-1]+nrecv[i-1];
  }
  // Fill up send_keys array with key-value pairs using the linked list
  int count;
  for (i=0; i<p_nprocs; i++) {
    int offset = s_offsets[i];
    if (ndest[i] > 0) {
      int count = offset;
      int j = ltop[i];
      while (j >= 0) {
        send_keys[count] = keys[j];
        j = ldest[j];
        count++;
      }
    }
  }
  // Send buffer is full so distribute contents to all processors
  ierr = MPI_Alltoallv(&send_keys[0], &ndest[0], &s_offsets[0], MPI_INT,
      &recv_keys[0], &nrecv[0], &r_offsets[0], MPI_INT, p_comm);
  // keys are available, so find corresponding values. Need to pass through data
  // twice, once to evaluate how many values are being returned and once to set
  // up return data structures
  std::multimap<int, int>::iterator it;
  int lo, hi;
  for (i=0; i<p_nprocs; i++) {
    lo = r_offsets[i];
    hi = lo+nrecv[i];
    ndest[i] = 0;
    for (j=lo; j<hi; j++) {
      it = p_umap.find(recv_keys[j]);
      if (it != p_umap.end()) {
        while (it != p_umap.upper_bound(recv_keys[j])) {
          ndest[i]++;
          it++;
        }
      }
    }
  }
  s_offsets[0] = 0;
  size = ndest[0];
  for (i=1; i<p_nprocs; i++) {
    size += ndest[i];
    s_offsets[i] = s_offsets[i-1] + ndest[i-1];
  }
  for (i=0; i<p_nprocs; i++) {
    ndest[i] = 2*ndest[i];
    s_offsets[i] = 2*s_offsets[i];
  }
  std::vector<int> send_values(2*size);
  // Pack returning data into send buffer
  count = 0;
  for (i=0; i<p_nprocs; i++) {
    lo = r_offsets[i];
    hi = lo+nrecv[i];
    for (j=lo; j<hi; j++) {
      it = p_umap.find(recv_keys[j]);
      if (it != p_umap.end()) {
        while (it != p_umap.upper_bound(recv_keys[j])) {
          send_values[2*count] = recv_keys[j];
          send_values[2*count+1] = it->second;
          count++;
          it++;
        }
      }
    }
  }
  // Send dimensions of returning blocks to all processors
  ierr = MPI_Alltoall(&ndest[0], one, MPI_INT, &nrecv[0], one, MPI_INT, p_comm);

  // Evaluate offsets for returning data
  r_offsets[0] = 0;
  rsize = nrecv[0];
  for (i=1; i<p_nprocs; i++) {
    r_offsets[i] = r_offsets[i-1] + nrecv[i-1];
    rsize += nrecv[i];
  }
  std::vector<int> ret_values(rsize);
  // send values back to requesting processor using all-to-all
  ierr = MPI_Alltoallv(&send_values[0], &ndest[0], &s_offsets[0], MPI_INT,
      &ret_values[0], &nrecv[0], &r_offsets[0], MPI_INT, p_comm);
  ierr = MPI_Barrier(p_comm);

  // repack data into output
  keys.clear();
  values.clear();
  size = rsize/2;
  for (i=0; i<size; i++) {
    keys.push_back(ret_values[2*i]);
    values.push_back(ret_values[2*i+1]);
  }
#else
  // Create a global array to count how many values are coming from each
  // processor and use this to create a set of offsets
  int g_offset = GA_Create_handle();
  int dims = p_nprocs;
  int one = 1;
  int blocks = 1;
  GA_Set_data(g_offset,one,&dims,C_INT);
  GA_Set_chunk(g_offset, &blocks);
  GA_Set_pgroup(g_offset, p_GAgrp);
  GA_Allocate(g_offset);
  GA_Zero(g_offset);

  // Get offsets on remote processors
  std::vector<int> r_offset(p_nprocs);
  for (j=0; j<p_nprocs; j++) {
    i = (j+p_me)%p_nprocs;
    if (ndest[i] > 0) {
      r_offset[i] = NGA_Read_inc(g_offset, &i, ndest[i]);
    } else {
      r_offset[i] = 0;
    }
    r_offset[i] = 2*r_offset[i];
  }
  GA_Pgroup_sync(p_GAgrp);
  // Copy values from global array to local array and evaluate global offsets
  std::vector<int> numValues(p_nprocs);
  int lo, hi;
  lo = 0;
  hi = p_nprocs-1;
  if (p_me == 0) {
    NGA_Get(g_offset,&lo,&hi,&numValues[0],&one);
  } else {
    for (i=0; i<p_nprocs; i++) {
      numValues[i] = 0;
    }
  }
  char plus[2];
  strcpy(plus,"+");
  GA_Pgroup_igop(p_GAgrp, &numValues[0], p_nprocs, plus);
  int totalVals = 0;
  for (i=0; i<p_nprocs; i++) {
    r_offset[i] += totalVals;
    totalVals += 2*numValues[i];
  }
  if (totalVals == 0) return;

  // Create a global array that can hold all key values
  int g_data = GA_Create_handle();
  dims = totalVals+1;
  GA_Set_data(g_data,one,&dims,C_INT);
  GA_Set_pgroup(g_data,p_GAgrp);
  std::vector<int> mapc(p_nprocs);
  mapc[0] = 0;
  for (i=1; i<p_nprocs; i++) {
    mapc[i] = mapc[i-1] + 2*numValues[i-1];
  }
  blocks = p_nprocs;
  GA_Set_irreg_distr(g_data,&mapc[0],&blocks);
  GA_Allocate(g_data);

  // Repack keys and send them to the processor that owns the corresponding
  // values
  int *data_keys;
  int ncnt;
  for (i=0; i<p_nprocs; i++) {
    j = ltop[i];
    ncnt = 0;
    if (j >= 0) {
      data_keys = new int[2*ndest[i]];
      while (j >= 0) {
        data_keys[2*ncnt] = keys[j];
        data_keys[2*ncnt+1] = p_me;
        j = ldest[j];
        ncnt++;
      }
      lo = r_offset[i];
      hi = lo + 2*ndest[i] - 1;
      NGA_Put(g_data, &lo, &hi, data_keys, &one);
      delete [] data_keys;
    }
  }
  GA_Pgroup_sync(p_GAgrp);

  // Keys are now all on processor holding the values, along with information on
  // the processor requesting the values.
  int nval = numValues[p_me];
  lo = mapc[p_me];
  hi = lo + 2*numValues[p_me] - 1;
  int ld;
  std::vector<int> data_pairs;
  if (lo<=hi) {
    NGA_Access(g_data,&lo,&hi,&data_keys,&ld);
    std::multimap<int,int>::iterator it;
    for (i=0; i<nval; i++) {
      it = p_umap.find(data_keys[2*i]);
      if (it != p_umap.end()) {
        while(it != p_umap.upper_bound(data_keys[2*i])) {
          data_pairs.push_back(data_keys[2*i]);
          data_pairs.push_back(data_keys[2*i+1]);
          data_pairs.push_back(it->second);
          it++;
        }
      } else {
        printf("p[%d] (index_hash) key not found: %d\n",p_me,data_keys[2*i]);
      }
    }
    NGA_Release(g_data,&lo,&hi);
  }
  GA_Destroy(g_data);

  // Now need to repack data and send it back to requesting processor
  size = data_pairs.size()/3;
  std::vector<int> lreturn(size);
  // initialize all arrays
  for (i=0; i<p_nprocs; i++) {
    ndest[i] = 0;
    ltop[i] = -1;
  }
  for (i=0; i<size; i++) {
    lreturn[i] = -1;
  }
  // create linked list of return values
  for (i=0; i<size; i++) {
    j = data_pairs[3*i+1];
    ndest[j]++;
    lreturn[i] = ltop[j];
    ltop[j] = i;
  }

  // Evaluate offsets for returning values
  GA_Zero(g_offset);

  // Get offsets on remote processors
  for (j=0; j<p_nprocs; j++) {
    i = (j+p_me)%p_nprocs;
    if (ndest[i] > 0) {
      r_offset[i] = NGA_Read_inc(g_offset, &i, ndest[i]);
    } else {
      r_offset[i] = 0;
    }
    r_offset[i] = 2*r_offset[i];
  }
  GA_Pgroup_sync(p_GAgrp);

  // Copy values from global array to local array and evaluate global offsets
  lo = 0;
  hi = p_nprocs-1;
  if (p_me == 0) {
    NGA_Get(g_offset,&lo,&hi,&numValues[0],&one);
  } else {
    for (i=0; i<p_nprocs; i++) {
      numValues[i] = 0;
    }
  }
  GA_Pgroup_igop(p_GAgrp, &numValues[0], p_nprocs, plus);
  GA_Destroy(g_offset);
  totalVals = 0;
  for (i=0; i<p_nprocs; i++) {
    r_offset[i] += totalVals;
    totalVals += 2*numValues[i];
  }

  // Create a global array that can hold all key-value pairs
  g_data = GA_Create_handle();
  dims = totalVals+1;
  GA_Set_data(g_data,one,&dims,C_INT);
  GA_Set_pgroup(g_data,p_GAgrp);
  mapc[0] = 0;
  for (i=1; i<p_nprocs; i++) {
    mapc[i] = mapc[i-1] + 2*numValues[i-1];
  }
  blocks = p_nprocs;
  GA_Set_irreg_distr(g_data,&mapc[0],&blocks);
  GA_Allocate(g_data);

  // Repack key-value pairs and send them to the processor that requested the
  // values
  for (i=0; i<p_nprocs; i++) {
    j = ltop[i];
    ncnt = 0;
    if (j >= 0) {
      data_keys = new int[2*ndest[i]];
      while (j >= 0) {
        data_keys[2*ncnt] = data_pairs[3*j];
        data_keys[2*ncnt+1] = data_pairs[3*j+2];
        j = lreturn[j];
        ncnt++;
      }
      lo = r_offset[i];
      hi = lo + 2*ndest[i] - 1;
      NGA_Put(g_data, &lo, &hi, data_keys, &one);
      delete [] data_keys;
    }
  }
  GA_Pgroup_sync(p_GAgrp);

  // Data has been returned to requesting process. Unpack it and fill the keys
  // and values arrays
  nval = numValues[p_me];
  lo = mapc[p_me];
  hi = lo + 2*numValues[p_me] - 1;
  keys.clear();
  values.clear();
  if (lo<=hi) {
    NGA_Access(g_data,&lo,&hi,&data_keys,&ld);
    for (i=0; i<nval; i++) {
      keys.push_back(data_keys[2*i]);
      values.push_back(data_keys[2*i+1]);
    }
    NGA_Release(g_data,&lo,&hi);
  }
  GA_Destroy(g_data);
#endif
}

// get values corresponding to a list of keys from the hash map where key is a
// pair of integers
// @param keys list of integer-pair keys
// @param values returned list of values corresponding to the list of keys
void GlobalIndexHashMap::getValues(std::vector<std::pair<int,int> > &keys,
    std::vector<int> &values)
{
  // Need to distribute keys to processors that hold the corresponding values.
  // Start by constructing a linked list of where each key needs to go
  int i, j;
  int size = keys.size();
  // ndest[idx] number of values to send to processor idx, nrecv[idx] is number
  // of values received from processor idx
  std::vector<int> ndest(p_nprocs);
  std::vector<int> nrecv(p_nprocs);
  // arrays for linked list
  std::vector<int> ltop(p_nprocs);
  std::vector<int> ldest(size);
  // initialize all arrays
  for (i=0; i<p_nprocs; i++) {
    ndest[i] = 0;
    nrecv[i] = 0;
    ltop[i] = -1;
  }
  for (i=0; i<size; i++) {
    ldest[i] = -1;
  }
  // create linked list
  int hash;
  for (i=0; i<size; i++) {
    hash = pairHashValue(keys[i]);
    ndest[hash]++;
    ldest[i] = ltop[hash];
    ltop[hash] = i;
  }
#ifdef HASH_WITH_MPI
  // send keys to processors based on linked list. Start by evaluating how much
  // data will be received from other processors using an all-to-all call
  int ierr;
  int one = 1;
  ierr = MPI_Alltoall(&ndest[0], one, MPI_INT, &nrecv[0], one, MPI_INT, p_comm);
  // nrecv now contains the number of keys that will be received from other
  // processors. Use this information to set up send and receive along with
  // their offsets for an all-to-all-v call
  int rsize = 0;
  for (i=0; i<p_nprocs; i++) {
    rsize += nrecv[i];
  }
  std::vector<int> send_keys(2*size);
  std::vector<int> recv_keys(2*rsize);
  std::vector<int> s_offsets(p_nprocs);
  std::vector<int> r_offsets(p_nprocs);
  s_offsets[0] = 0;
  r_offsets[0] = 0;
  for (i=1; i<p_nprocs; i++) {
    s_offsets[i] = s_offsets[i-1]+2*ndest[i-1];
    r_offsets[i] = r_offsets[i-1]+2*nrecv[i-1];
  }
  for (i=0; i<p_nprocs; i++) {
    ndest[i] = 2*ndest[i];
    nrecv[i] = 2*nrecv[i];
  }
  // Fill up send_keys array with key-value pairs using the linked list
  int count;
  for (i=0; i<p_nprocs; i++) {
    int offset = s_offsets[i];
    if (ndest[i] > 0) {
      int count = offset;
      int j = ltop[i];
      while (j >= 0) {
        send_keys[count] = keys[j].first;
        send_keys[count+1] = keys[j].second;
        j = ldest[j];
        count += 2;
      }
    }
  }
  // Send buffer is full so distribute contents to all processors
  ierr = MPI_Alltoallv(&send_keys[0], &ndest[0], &s_offsets[0], MPI_INT,
      &recv_keys[0], &nrecv[0], &r_offsets[0], MPI_INT, p_comm);
  // keys are available, so find corresponding values. Need to pass through data
  // twice, once to evaluate how many values are being returned and once to set
  // up return data structures
  std::multimap<std::pair<int,int>,int>::iterator it;
  std::pair<int,int> key;
  int lo, hi;
  for (i=0; i<p_nprocs; i++) {
    lo = r_offsets[i]/2;
    hi = lo+nrecv[i]/2;
    ndest[i] = 0;
    for (j=lo; j<hi; j++) {
      key = std::pair<int,int>(recv_keys[2*j],recv_keys[2*j+1]);
      it = p_pmap.find(key);
      if (it != p_pmap.end()) {
        while (it != p_pmap.upper_bound(key)) {
          ndest[i]++;
          it++;
        }
      }
    }
  }
  s_offsets[0] = 0;
  size = ndest[0];
  for (i=1; i<p_nprocs; i++) {
    size += ndest[i];
    s_offsets[i] = s_offsets[i-1] + ndest[i-1];
  }
  for (i=0; i<p_nprocs; i++) {
    ndest[i] = 3*ndest[i];
    s_offsets[i] = 3*s_offsets[i];
  }
  std::vector<int> send_values(3*size);
  // Pack returning data into send buffer
  count = 0;
  for (i=0; i<p_nprocs; i++) {
    lo = r_offsets[i]/2;
    hi = lo+nrecv[i]/2;
    for (j=lo; j<hi; j++) {
      key = std::pair<int,int>(recv_keys[2*j],recv_keys[2*j+1]);
      it = p_pmap.find(key);
      if (it != p_pmap.end()) {
        while (it != p_pmap.upper_bound(key)) {
          send_values[3*count] = key.first;
          send_values[3*count+1] = key.second;
          send_values[3*count+2] = it->second;
          count++;
          it++;
        }
      }
    }
  }
  // Send dimensions of returning blocks to all processors
  ierr = MPI_Alltoall(&ndest[0], one, MPI_INT, &nrecv[0], one, MPI_INT, p_comm);

  // Evaluate offsets for returning data
  r_offsets[0] = 0;
  rsize = nrecv[0];
  for (i=1; i<p_nprocs; i++) {
    r_offsets[i] = r_offsets[i-1] + nrecv[i-1];
    rsize += nrecv[i];
  }
  std::vector<int> ret_values(rsize);
  // send values back to requesting processor using all-to-all
  ierr = MPI_Alltoallv(&send_values[0], &ndest[0], &s_offsets[0], MPI_INT,
      &ret_values[0], &nrecv[0], &r_offsets[0], MPI_INT, p_comm);
  ierr = MPI_Barrier(p_comm);

  // repack data into output
  keys.clear();
  values.clear();
  size = rsize/3;
  for (i=0; i<size; i++) {
    key = std::pair<int,int>(ret_values[3*i],ret_values[3*i+1]);
    keys.push_back(key);
    values.push_back(ret_values[3*i+2]);
  }
#else
  // Create a global array to count how many values are coming from each
  // processor and use this to create a set of offsets
  int g_offset = GA_Create_handle();
  int dims = p_nprocs;
  int one = 1;
  int blocks = 1;
  GA_Set_data(g_offset,one,&dims,C_INT);
  GA_Set_chunk(g_offset, &blocks);
  GA_Set_pgroup(g_offset, p_GAgrp);
  GA_Allocate(g_offset);
  GA_Zero(g_offset);

  // Get offsets on remote processors
  std::vector<int> r_offset(p_nprocs);
  for (j=0; j<p_nprocs; j++) {
    i = (j+p_me)%p_nprocs;
    if (ndest[i] > 0) {
      r_offset[i] = NGA_Read_inc(g_offset, &i, ndest[i]);
    } else {
      r_offset[i] = 0;
    }
    r_offset[i] = 3*r_offset[i];
  }
  GA_Pgroup_sync(p_GAgrp);
  // Copy values from global array to local array and evaluate global offsets
  std::vector<int> numValues(p_nprocs);
  int lo, hi;
  lo = 0;
  hi = p_nprocs-1;
  if (p_me == 0) {
    NGA_Get(g_offset,&lo,&hi,&numValues[0],&one);
  } else {
    for (i=0; i<p_nprocs; i++) {
      numValues[i] = 0;
    }
  }
  char plus[2];
  strcpy(plus,"+");
  GA_Pgroup_igop(p_GAgrp, &numValues[0], p_nprocs, plus);
  int totalVals = 0;
  for (i=0; i<p_nprocs; i++) {
    r_offset[i] += totalVals;
    totalVals += 3*numValues[i];
  }
  if (totalVals == 0) return;

  // Create a global array that can hold all key values
  int g_data = GA_Create_handle();
  dims = totalVals+1;
  GA_Set_data(g_data,one,&dims,C_INT);
  GA_Set_pgroup(g_data,p_GAgrp);
  std::vector<int> mapc(p_nprocs);
  mapc[0] = 0;
  for (i=1; i<p_nprocs; i++) {
    mapc[i] = mapc[i-1] + 3*numValues[i-1];
  }
  blocks = p_nprocs;
  GA_Set_irreg_distr(g_data,&mapc[0],&blocks);
  GA_Allocate(g_data);

  // Repack keys and send them to the processor that owns the corresponding
  // values
  int *data_keys;
  int ncnt;
  for (i=0; i<p_nprocs; i++) {
    j = ltop[i];
    ncnt = 0;
    if (j >= 0) {
      data_keys = new int[3*ndest[i]];
      while (j >= 0) {
        data_keys[3*ncnt] = keys[j].first;
        data_keys[3*ncnt+1] = keys[j].second;
        data_keys[3*ncnt+2] = p_me;
        j = ldest[j];
        ncnt++;
      }
      lo = r_offset[i];
      hi = lo + 3*ndest[i] - 1;
      NGA_Put(g_data, &lo, &hi, data_keys, &one);
      delete [] data_keys;
    }
  }
  GA_Pgroup_sync(p_GAgrp);

  // Keys are now all on processor holding the values, along with information on
  // the processor requesting the values.
  int nval = numValues[p_me];
  lo = mapc[p_me];
  hi = lo + 3*numValues[p_me] - 1;
  int ld;
  std::vector<int> data_pairs;
  if (lo<=hi) {
    NGA_Access(g_data,&lo,&hi,&data_keys,&ld);
    std::multimap<std::pair<int,int>,int>::iterator it;
    for (i=0; i<nval; i++) {
      std::pair<int,int> key
        = std::pair<int,int>(data_keys[3*i],data_keys[3*i+1]);
      it = p_pmap.find(key);
      if (it != p_pmap.end()) {
        while(it != p_pmap.upper_bound(key)) {
          data_pairs.push_back(data_keys[3*i]);
          data_pairs.push_back(data_keys[3*i+1]);
          data_pairs.push_back(data_keys[3*i+2]);
          data_pairs.push_back(it->second);
          it++;
        }
      } else {
        printf("p[%d] (index_hash) key not found: < %d, %d>\n",
            p_me,data_keys[3*i],data_keys[3*i+1]);
      }
    }
    NGA_Release(g_data,&lo,&hi);
  }
  GA_Destroy(g_data);

  // Now need to repack data and send it back to requesting processor
  size = data_pairs.size()/4;
  std::vector<int> lreturn(size);
  // initialize all arrays
  for (i=0; i<p_nprocs; i++) {
    ndest[i] = 0;
    ltop[i] = -1;
  }
  for (i=0; i<size; i++) {
    lreturn[i] = -1;
  }
  // create linked list of return values
  for (i=0; i<size; i++) {
    j = data_pairs[4*i+2];
    ndest[j]++;
    lreturn[i] = ltop[j];
    ltop[j] = i;
  }

  // Evaluate offsets for returning values
  GA_Zero(g_offset);

  // Get offsets on remote processors
  for (j=0; j<p_nprocs; j++) {
    i = (j+p_me)%p_nprocs;
    if (ndest[i] > 0) {
      r_offset[i] = NGA_Read_inc(g_offset, &i, ndest[i]);
    } else {
      r_offset[i] = 0;
    }
    r_offset[i] = 3*r_offset[i];
  }
  GA_Pgroup_sync(p_GAgrp);

  // Copy values from global array to local array and evaluate global offsets
  lo = 0;
  hi = p_nprocs-1;
  if (p_me == 0) {
    NGA_Get(g_offset,&lo,&hi,&numValues[0],&one);
  } else {
    for (i=0; i<p_nprocs; i++) {
      numValues[i] = 0;
    }
  }
  GA_Pgroup_igop(p_GAgrp, &numValues[0], p_nprocs, plus);
  GA_Destroy(g_offset);
  totalVals = 0;
  for (i=0; i<p_nprocs; i++) {
    r_offset[i] += totalVals;
    totalVals += 3*numValues[i];
  }

  // Create a global array that can hold all key-value pairs
  g_data = GA_Create_handle();
  dims = totalVals+1;
  GA_Set_data(g_data,one,&dims,C_INT);
  GA_Set_pgroup(g_data,p_GAgrp);
  mapc[0] = 0;
  for (i=1; i<p_nprocs; i++) {
    mapc[i] = mapc[i-1] + 3*numValues[i-1];
  }
  blocks = p_nprocs;
  GA_Set_irreg_distr(g_data,&mapc[0],&blocks);
  GA_Allocate(g_data);

  // Repack key-value pairs and send them to the processor that requested the
  // values
  for (i=0; i<p_nprocs; i++) {
    j = ltop[i];
    ncnt = 0;
    if (j >= 0) {
      data_keys = new int[3*ndest[i]];
      while (j >= 0) {
        data_keys[3*ncnt] = data_pairs[4*j];
        data_keys[3*ncnt+1] = data_pairs[4*j+1];
        data_keys[3*ncnt+2] = data_pairs[4*j+3];
        j = lreturn[j];
        ncnt++;
      }
      lo = r_offset[i];
      hi = lo + 3*ndest[i] - 1;
      NGA_Put(g_data, &lo, &hi, data_keys, &one);
      delete [] data_keys;
    }
  }
  GA_Pgroup_sync(p_GAgrp);

  // Data has been returned to requesting process. Unpack it and fill the keys
  // and values arrays
  nval = numValues[p_me];
  lo = mapc[p_me];
  hi = lo + 3*numValues[p_me] - 1;
  keys.clear();
  values.clear();
  if (lo<=hi) {
    NGA_Access(g_data,&lo,&hi,&data_keys,&ld);
    for (i=0; i<nval; i++) {
      keys.push_back(std::pair<int,int>(data_keys[3*i],data_keys[3*i+1]));
      values.push_back(data_keys[3*i+2]);
    }
    NGA_Release(g_data,&lo,&hi);
  }
  GA_Destroy(g_data);
#endif
}

// hash function for indices. Maps the value of key into the interval [0,p_nprocs-1]
// where p_nprocs is the total number of processors
// @param key input value of key
// @return number between 0 and p_nprocs-1
int GlobalIndexHashMap::hashValue(int key)
{
  return key%p_nprocs;
}

// hash function for index pairs. Maps the value of key into the interval [0,p_nprocs-1]
// where p_nprocs is the total number of processors
// @param key input value of key
// @return number between 0 and nprocs-1
int GlobalIndexHashMap::pairHashValue(std::pair<int,int> key)
{
  int i1 = key.first;
  int i2 = key.second;
  // 1009 is prime number
  return ((i2*1009+i1)%p_nprocs);
}

} // hash_map
} // gridpack
