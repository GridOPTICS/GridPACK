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

#include <parallel/communicator.hpp>
#include "index_hash.hpp"

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
}

// Default destructor
GlobalIndexHashMap::~GlobalIndexHashMap(void)
{
}

// add key-value pairs to hash map
void GlobalIndexHashMap::addPairs(std::vector<std::pair<int,int> > &pairs)
{
  // Need to distribute key-value pairs between processors based on the value
  // returned by the hashValue function. Start by constructing a linked list of
  // where each pair needs to go
  int i;
  int size = pairs.size();
  // ndest[idx] number of values to send to processor idx, nrecv[idx] is number
  // of values received from processor idx
  int ndest[p_nprocs], nrecv[p_nprocs];
  // arrays for linked list
  int ltop[p_nprocs], ldest[size];
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
  // send data to processors based on linked list. Start by evaluating how much
  // data will be received from other processors using an all-to-all call
  int ierr;
  int one = 1;
  ierr = MPI_Alltoall(ndest, one, MPI_INT, nrecv, one, MPI_INT, p_comm);
  // nrecv now contains the number of pairs that will be received from other
  // processors. Use this information to set up send and receive along with
  // their offsets for a all-to-all-v call
  int rsize = 0;
  for (i=0; i<p_nprocs; i++) {
    rsize += nrecv[i];
    ndest[i] = 2*ndest[i];
    nrecv[i] = 2*nrecv[i];
  }
  int send_pair[2*size];
  int recv_pair[2*rsize];
  int s_offsets[p_nprocs];
  int r_offsets[p_nprocs];
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
  ierr = MPI_Alltoallv(send_pair, ndest, s_offsets, MPI_INT,
      recv_pair, nrecv, r_offsets, MPI_INT, p_comm);
  // key-value pairs are available, so set up local hash table
  p_umap.clear();
  int first, second;
  for (i=0; i<rsize; i++) {
    first = recv_pair[2*i];
    second = recv_pair[2*i+1];
    p_umap.insert(std::pair<int, int>(first,second));
  }
}

// get values corresponding to a list of keys from the hash map
void GlobalIndexHashMap::getValues(std::vector<int> &keys, std::vector<int> &values)
{
  // Need to distribute keys to processors that hold the corresponding values.
  // Start by constructing a linked list of where each key needs to go
  int i;
  int size = keys.size();
  // ndest[idx] number of values to send to processor idx, nrecv[idx] is number
  // of values received from processor idx
  int ndest[p_nprocs], nrecv[p_nprocs];
  // arrays for linked list
  int ltop[p_nprocs], ldest[size];
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
  // send keys to processors based on linked list. Start by evaluating how much
  // data will be received from other processors using an all-to-all call
  int ierr;
  int one = 1;
  ierr = MPI_Alltoall(ndest, one, MPI_INT, nrecv, one, MPI_INT, p_comm);
  // nrecv now contains the number of keys that will be received from other
  // processors. Use this information to set up send and receive along with
  // their offsets for a all-to-all-v call
  int rsize = 0;
  for (i=0; i<p_nprocs; i++) {
    rsize += nrecv[i];
    ndest[i] = ndest[i];
    nrecv[i] = nrecv[i];
  }
  int send_keys[size];
  int recv_keys[rsize];
  int s_offsets[p_nprocs];
  int r_offsets[p_nprocs];
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
  ierr = MPI_Alltoallv(send_keys, ndest, s_offsets, MPI_INT,
      recv_keys, nrecv, r_offsets, MPI_INT, p_comm);
  // keys are available, so find corresponding values
  boost::unordered_map<int, int>::iterator it;
  int send_values[rsize];
  for (i=0; i<rsize; i++) {
    it = p_umap.find(recv_keys[i]);
    if (it != p_umap.end()) {
      send_values[i] = it->second;
    } else {
      // TODO: Some kind of error
    }
  }
  // send values back to requesting processor using all-to-all
  int recv_values[size];
  ierr = MPI_Alltoallv(send_values, nrecv, r_offsets, MPI_INT,
      recv_values, ndest, s_offsets, MPI_INT, p_comm);
  // values in recv_values buffer into values vector
  values.clear();
  values.reserve(size);
  for (i=0; i<p_nprocs; i++) {
    int offset = s_offsets[i];
    if (ndest[i] > 0) {
      int count = offset;
      int j = ltop[i];
      while (j >= 0) {
        values[j] = recv_values[count];
        j = ldest[j];
        count++;
      }
    }
  }
}

// hash function for indices. Maps the value of key into the interval [0,p_nprocs-1]
// where p_nprocs is the total number of processors
int GlobalIndexHashMap::hashValue(int key)
{
  return key%p_nprocs;
}

} // hash_map
} // gridpack
